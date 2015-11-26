// provide render functionality to the example. The code is based on the
// D3D tutorials, and is basically just a quick hack to get things
// working. Not really nicely written...

#define STRICT
#define D3D_OVERLOADS
#include <math.h>
#include <windows.h>
#include <stdio.h>
#include <d3d.h>
#include "resource.h"
#include "geo.h"

static LPDIRECTDRAW         g_pDD1           = NULL;
static LPDIRECTDRAW4        g_pDD4           = NULL;
static LPDIRECTDRAWSURFACE4 g_pddsPrimary    = NULL;
static LPDIRECTDRAWSURFACE4 g_pddsBackBuffer = NULL;
static LPDIRECTDRAWSURFACE4 g_pddsZBuffer    = NULL;
static LPDIRECT3D3          g_pD3D           = NULL;
static LPDIRECT3DDEVICE3    g_pd3dDevice     = NULL;
static LPDIRECT3DVIEWPORT3  g_pvViewport     = NULL;
static RECT                 g_rcScreenRect;
static RECT                 g_rcViewportRect;

static HWND MainWin;

static BOOL g_bActive  = FALSE;
static BOOL g_bReady   = FALSE;

LRESULT CALLBACK WndProc( HWND, UINT, WPARAM, LPARAM );
HRESULT CreateEverything( HWND );
HRESULT Initialize3DEnvironment( HWND, GUID*, const GUID* );
HRESULT Cleanup3DEnvironment();
HRESULT Render3DEnvironment();
VOID    OnMove( INT, INT );
HRESULT ShowFrame();
HRESULT RestoreSurfaces();

HRESULT App_FrameMove( LPDIRECT3DDEVICE3, FLOAT );

LPDIRECT3DMATERIAL3 g_pmtrlObjectMtrl = NULL;
D3DVERTEX           g_pvTriangleVertices[6];


DWORD VertexType=D3DFVF_XYZ | D3DFVF_DIFFUSE;
typedef struct {
	D3DVALUE x,y,z;
	DWORD diffcolor;
} MyVertex;
MyVertex Vertex[24];


HRESULT App_InitDeviceObjects( LPDIRECT3DDEVICE3 pd3dDevice,
							  LPDIRECT3DVIEWPORT3 pvViewport )
{
    int i;
    for (i=0;i<8;i++)   Vertex[i].diffcolor=D3DRGB(1,0,0);
    for (i=8;i<16;i++)  Vertex[i].diffcolor=D3DRGB(0,1,0);
    for (i=16;i<24;i++) Vertex[i].diffcolor=D3DRGB(0,0,1);
	
    LPDIRECT3D3 pD3D;
    if( FAILED( pd3dDevice->GetDirect3D( &pD3D ) ) ) return E_FAIL;
    pD3D->Release();
	
    if( FAILED( pD3D->CreateMaterial( &g_pmtrlObjectMtrl, NULL ) ) ) return E_FAIL;
	
    D3DMATERIAL       mtrl;
    D3DMATERIALHANDLE hmtrl;
    ZeroMemory( &mtrl, sizeof(D3DMATERIAL) );
    mtrl.dwSize       = sizeof(D3DMATERIAL);
    mtrl.dcvAmbient.r = 0.0f;
    mtrl.dcvAmbient.g = 0.0f;
    mtrl.dcvAmbient.b = 0.0f;
    g_pmtrlObjectMtrl->SetMaterial( &mtrl );
	
    g_pmtrlObjectMtrl->GetHandle( pd3dDevice, &hmtrl );
    pd3dDevice->SetLightState(  D3DLIGHTSTATE_MATERIAL, hmtrl );
	
    pd3dDevice->SetLightState(  D3DLIGHTSTATE_AMBIENT,  0xffffffff );
	
    D3DMATRIX mat;
	mat._11 = mat._22 = mat._33 = mat._44 = 1.0f;
	mat._12 = mat._13 = mat._14 = mat._41 = 0.0f;
	mat._21 = mat._23 = mat._24 = mat._42 = 0.0f;
	mat._31 = mat._32 = mat._34 = mat._43 = 0.0f;
	
    D3DMATRIX matWorld = mat;
    pd3dDevice->SetTransform( D3DTRANSFORMSTATE_WORLD, &matWorld );
	
    D3DMATRIX matView = mat;
    matView._43 = 10.0f;
    pd3dDevice->SetTransform( D3DTRANSFORMSTATE_VIEW, &matView );
	
    D3DMATRIX matProj = mat;
    matProj._11 =  2.0f;
    matProj._22 =  2.0f;
    matProj._34 =  1.0f;
    matProj._43 = -1.0f;
    matProj._44 =  0.0f;
    pd3dDevice->SetTransform( D3DTRANSFORMSTATE_PROJECTION, &matProj );
	
    return S_OK;
}

VOID App_DeleteDeviceObjects( LPDIRECT3DDEVICE3 pd3dDevice, 
							 LPDIRECT3DVIEWPORT3 pvViewport )
{
    if( g_pmtrlObjectMtrl ) g_pmtrlObjectMtrl->Release();
	g_pmtrlObjectMtrl = NULL;
}

HINSTANCE hInst;
void main();
MSG  msg;

INT WINAPI WinMain( HINSTANCE _hInst, HINSTANCE, LPSTR strCmdLine, INT )
{
	hInst=_hInst;
	main();
	return msg.wParam;
}



LRESULT CALLBACK WndProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam )
{
    switch( uMsg )
    {
	case WM_PAINT:
		if (hWnd==MainWin) ShowFrame();
		break;
		
	case WM_MOVE:
		if( g_bActive && g_bReady )
			OnMove( (SHORT)LOWORD(lParam), (SHORT)HIWORD(lParam) );
		break;
		
	case WM_SIZE:
		if( SIZE_MAXHIDE==wParam || SIZE_MINIMIZED==wParam )
			g_bActive = FALSE;
		else g_bActive = TRUE;
		
		if( g_bActive && g_bReady )
		{
			g_bReady = FALSE;
			if( FAILED( CreateEverything( MainWin ) ) )
				DestroyWindow( MainWin );
			g_bReady = TRUE;
		}
		break;
		
	case WM_GETMINMAXINFO:
		((MINMAXINFO*)lParam)->ptMinTrackSize.x = 100;
		((MINMAXINFO*)lParam)->ptMinTrackSize.y = 100;
		break;
		
	case WM_CLOSE:
		DestroyWindow( hWnd );
		return 0;
        
	case WM_DESTROY:
		Cleanup3DEnvironment();
		PostQuitMessage(0);
		return 0L;
    }
	
    return DefWindowProc( hWnd, uMsg, wParam, lParam );
}




HRESULT CreateEverything( HWND hWnd )
{
	if( FAILED( Cleanup3DEnvironment() ) )
		return E_FAIL;
				
	if( SUCCEEDED( Initialize3DEnvironment( hWnd, NULL,	&IID_IDirect3DHALDevice ) ) ) return S_OK;
	
	Cleanup3DEnvironment();
	
	if( SUCCEEDED( Initialize3DEnvironment( hWnd, NULL,	&IID_IDirect3DRGBDevice ) ) ) return S_OK;
	
	return E_FAIL;
}



static HRESULT WINAPI EnumZBufferCallback( DDPIXELFORMAT* pddpf,
										  VOID* pddpfDesired ){
    if( pddpf->dwFlags == DDPF_ZBUFFER )    {
        memcpy( pddpfDesired, pddpf, sizeof(DDPIXELFORMAT) ); 
        return D3DENUMRET_CANCEL;
    } 
    return D3DENUMRET_OK;
} 



HRESULT Initialize3DEnvironment( HWND hWnd, GUID* pDriverGUID, const GUID* pDeviceGUID )
{
	HRESULT hr;
	
	hr = DirectDrawCreate( pDriverGUID, &g_pDD1, NULL );
	if( FAILED( hr ) ) return hr;
	
	hr = g_pDD1->QueryInterface( IID_IDirectDraw4, (VOID**)&g_pDD4 );
	if( FAILED( hr ) ) return hr;
	
    hr = g_pDD4->SetCooperativeLevel( hWnd, DDSCL_NORMAL );
	if( FAILED( hr ) ) return hr;
	
	DDSURFACEDESC2 ddsd;
	ZeroMemory( &ddsd, sizeof(DDSURFACEDESC2) );
	ddsd.dwSize         = sizeof(DDSURFACEDESC2);
	ddsd.dwFlags        = DDSD_CAPS;
	ddsd.ddsCaps.dwCaps = DDSCAPS_PRIMARYSURFACE;
	
	hr = g_pDD4->CreateSurface( &ddsd, &g_pddsPrimary, NULL );
	if( FAILED( hr ) ) return hr;
	
	ddsd.dwFlags        = DDSD_WIDTH | DDSD_HEIGHT | DDSD_CAPS;
	ddsd.ddsCaps.dwCaps = DDSCAPS_OFFSCREENPLAIN | DDSCAPS_3DDEVICE;
	
	GetClientRect( hWnd, &g_rcScreenRect );
	GetClientRect( hWnd, &g_rcViewportRect );
	ClientToScreen( hWnd, (POINT*)&g_rcScreenRect.left );
	ClientToScreen( hWnd, (POINT*)&g_rcScreenRect.right );
	ddsd.dwWidth  = g_rcScreenRect.right - g_rcScreenRect.left;
	ddsd.dwHeight = g_rcScreenRect.bottom - g_rcScreenRect.top;
	
	hr = g_pDD4->CreateSurface( &ddsd, &g_pddsBackBuffer, NULL );
	if( FAILED( hr ) ) return hr;
	
	LPDIRECTDRAWCLIPPER pcClipper;
	hr = g_pDD4->CreateClipper( 0, &pcClipper, NULL );
	if( FAILED( hr ) ) return hr;
	
	pcClipper->SetHWnd( 0, hWnd );
	g_pddsPrimary->SetClipper( pcClipper );
	pcClipper->Release();
	
    g_pDD4->QueryInterface( IID_IDirect3D3, (VOID**)&g_pD3D );
    if( FAILED( hr) ) return hr;
	
    DDPIXELFORMAT ddpfZBuffer;
    g_pD3D->EnumZBufferFormats( *pDeviceGUID, 
		EnumZBufferCallback, (VOID*)&ddpfZBuffer );
	
    if( sizeof(DDPIXELFORMAT) != ddpfZBuffer.dwSize ) return E_FAIL;
	
    ddsd.dwFlags        = DDSD_CAPS|DDSD_WIDTH|DDSD_HEIGHT|DDSD_PIXELFORMAT;
    ddsd.ddsCaps.dwCaps = DDSCAPS_ZBUFFER;
    ddsd.dwWidth        = g_rcScreenRect.right - g_rcScreenRect.left;
    ddsd.dwHeight       = g_rcScreenRect.bottom - g_rcScreenRect.top;
    memcpy( &ddsd.ddpfPixelFormat, &ddpfZBuffer, sizeof(DDPIXELFORMAT) );
	
    if( IsEqualIID( *pDeviceGUID, IID_IDirect3DHALDevice ) ) ddsd.ddsCaps.dwCaps |= DDSCAPS_VIDEOMEMORY;
    else                                                     ddsd.ddsCaps.dwCaps |= DDSCAPS_SYSTEMMEMORY;
	
    if( FAILED( hr = g_pDD4->CreateSurface( &ddsd, &g_pddsZBuffer, NULL ) ) ) return hr;
    if( FAILED( hr = g_pddsBackBuffer->AddAttachedSurface( g_pddsZBuffer ) ) ) return hr;
	
	ddsd.dwSize = sizeof(DDSURFACEDESC2);
	g_pDD4->GetDisplayMode( &ddsd );
	if( ddsd.ddpfPixelFormat.dwRGBBitCount <= 8 ) return DDERR_INVALIDMODE;
	
    hr = g_pD3D->CreateDevice( IID_IDirect3DHALDevice, g_pddsBackBuffer, &g_pd3dDevice, NULL );
	if( FAILED( hr ) )
	{
		hr = g_pD3D->CreateDevice( IID_IDirect3DRGBDevice, g_pddsBackBuffer,
			&g_pd3dDevice, NULL );
		if( FAILED( hr ) )
			return hr;
	}
	
    D3DVIEWPORT2 vdData;
    ZeroMemory( &vdData, sizeof(D3DVIEWPORT2) );
    vdData.dwSize       = sizeof(D3DVIEWPORT2);
	vdData.dwWidth      = g_rcScreenRect.right - g_rcScreenRect.left;
	vdData.dwHeight     = g_rcScreenRect.bottom - g_rcScreenRect.top;
    vdData.dvClipX      = -1.0f;
    vdData.dvClipWidth  = 2.0f;
    vdData.dvClipY      = 1.0f;
    vdData.dvClipHeight = 2.0f;
    vdData.dvMaxZ       = 1.0f;
	
    hr = g_pD3D->CreateViewport( &g_pvViewport, NULL );
	if( FAILED( hr ) )	return hr;
	
    g_pd3dDevice->AddViewport( g_pvViewport );
    g_pvViewport->SetViewport2( &vdData );
    g_pd3dDevice->SetCurrentViewport( g_pvViewport );
    g_pd3dDevice->SetRenderState( D3DRENDERSTATE_ZENABLE, TRUE );
	
	return App_InitDeviceObjects( g_pd3dDevice, g_pvViewport );
}

HRESULT Cleanup3DEnvironment()
{
	App_DeleteDeviceObjects( g_pd3dDevice, g_pvViewport );
	
    if( g_pvViewport )     g_pvViewport->Release();
	if( g_pD3D )           g_pD3D->Release();
	if( g_pddsBackBuffer ) g_pddsBackBuffer->Release();
	if( g_pddsPrimary )    g_pddsPrimary->Release();
	if( g_pDD4 )           g_pDD4->Release();
	
    if( g_pd3dDevice )
        if( 0 < g_pd3dDevice->Release() )
			return E_FAIL;
		
	if( g_pDD1 )
		if( 0 < g_pDD1->Release() ) return E_FAIL;
			
	g_pvViewport     = NULL;
	g_pd3dDevice     = NULL;
	g_pD3D           = NULL;
	g_pddsBackBuffer = NULL;
	g_pddsPrimary    = NULL;
	g_pDD4           = NULL;
	g_pDD1           = NULL;
			
	return S_OK;
}


HRESULT App_Render( LPDIRECT3DDEVICE3, LPDIRECT3DVIEWPORT3, D3DRECT*);

HRESULT Render3DEnvironment()
{
    App_Render( g_pd3dDevice, g_pvViewport, (D3DRECT*)&g_rcViewportRect );
    if( DDERR_SURFACELOST == ShowFrame() ) RestoreSurfaces();
    return S_OK;
}

HRESULT ShowFrame()
{
	if( NULL == g_pddsPrimary )	return E_FAIL;
    return g_pddsPrimary->Blt( &g_rcScreenRect, g_pddsBackBuffer, &g_rcViewportRect, DDBLT_WAIT, NULL );
}


HRESULT RestoreSurfaces()
{
    if( g_pddsPrimary )
        if( g_pddsPrimary->IsLost() )
            g_pddsPrimary->Restore();
		
	if( g_pddsBackBuffer )
		if( g_pddsBackBuffer->IsLost() )
			g_pddsBackBuffer->Restore();
			
	return S_OK;
}


VOID OnMove( INT x, INT y )
{
	DWORD dwWidth  = g_rcScreenRect.right - g_rcScreenRect.left;
	DWORD dwHeight = g_rcScreenRect.bottom - g_rcScreenRect.top;
    SetRect( &g_rcScreenRect, x, y, x + dwWidth, y + dwHeight );
}

HACCEL hAccel;
int InitRender(DL_dyna_system& system, MyCube* cube[], int ncubes){
    WNDCLASS wndClass = { CS_HREDRAW | CS_VREDRAW, WndProc, 0, 0, hInst,
		LoadIcon( hInst, MAKEINTRESOURCE(IDI_MAIN_ICON)),
		LoadCursor(NULL, IDC_ARROW), 
		(HBRUSH)GetStockObject(WHITE_BRUSH), NULL,
		TEXT("Render Window") };
    RegisterClass( &wndClass );
	
    MainWin = CreateWindow( TEXT("Render Window"),
		TEXT("Dynamo demo"),
		WS_OVERLAPPEDWINDOW,
		100, 10, 600, 600,
		0L, 0L, hInst, 0L );
    ShowWindow( MainWin, SW_SHOWNORMAL );
    UpdateWindow( MainWin );
	
    hAccel = LoadAccelerators( hInst, MAKEINTRESOURCE(IDR_MAIN_ACCEL) );
	
	if( FAILED( CreateEverything( MainWin ) ) )  return 0;
	return 1;
}

void SetVertices(DL_point& pos, DL_matrix& orient){
	Vertex[16].x=Vertex[12].x=Vertex[0].x=(float)(pos.x-orient.c0.x-orient.c1.x+orient.c2.x);
	Vertex[16].y=Vertex[12].y=Vertex[0].y=(float)(pos.y-orient.c0.y-orient.c1.y+orient.c2.y);
	Vertex[16].z=Vertex[12].z=Vertex[0].z=(float)(pos.z-orient.c0.z-orient.c1.z+orient.c2.z);
	
	Vertex[17].x=Vertex[8].x=Vertex[1].x=(float)(pos.x+orient.c0.x-orient.c1.x+orient.c2.x);
	Vertex[17].y=Vertex[8].y=Vertex[1].y=(float)(pos.y+orient.c0.y-orient.c1.y+orient.c2.y);
	Vertex[17].z=Vertex[8].z=Vertex[1].z=(float)(pos.z+orient.c0.z-orient.c1.z+orient.c2.z);
	
	Vertex[21].x=Vertex[9].x=Vertex[2].x=(float)(pos.x+orient.c0.x-orient.c1.x-orient.c2.x);
	Vertex[21].y=Vertex[9].y=Vertex[2].y=(float)(pos.y+orient.c0.y-orient.c1.y-orient.c2.y);
	Vertex[21].z=Vertex[9].z=Vertex[2].z=(float)(pos.z+orient.c0.z-orient.c1.z-orient.c2.z);
	
	Vertex[20].x=Vertex[13].x=Vertex[3].x=(float)(pos.x-orient.c0.x-orient.c1.x-orient.c2.x);
	Vertex[20].y=Vertex[13].y=Vertex[3].y=(float)(pos.y-orient.c0.y-orient.c1.y-orient.c2.y);
	Vertex[20].z=Vertex[13].z=Vertex[3].z=(float)(pos.z-orient.c0.z-orient.c1.z-orient.c2.z);
	
	Vertex[19].x=Vertex[15].x=Vertex[4].x=(float)(pos.x-orient.c0.x+orient.c1.x+orient.c2.x);
	Vertex[19].y=Vertex[15].y=Vertex[4].y=(float)(pos.y-orient.c0.y+orient.c1.y+orient.c2.y);
	Vertex[19].z=Vertex[15].z=Vertex[4].z=(float)(pos.z-orient.c0.z+orient.c1.z+orient.c2.z);
	
	Vertex[18].x=Vertex[11].x=Vertex[5].x=(float)(pos.x+orient.c0.x+orient.c1.x+orient.c2.x);
	Vertex[18].y=Vertex[11].y=Vertex[5].y=(float)(pos.y+orient.c0.y+orient.c1.y+orient.c2.y);
	Vertex[18].z=Vertex[11].z=Vertex[5].z=(float)(pos.z+orient.c0.z+orient.c1.z+orient.c2.z);
	
	Vertex[22].x=Vertex[10].x=Vertex[6].x=(float)(pos.x+orient.c0.x+orient.c1.x-orient.c2.x);
	Vertex[22].y=Vertex[10].y=Vertex[6].y=(float)(pos.y+orient.c0.y+orient.c1.y-orient.c2.y);
	Vertex[22].z=Vertex[10].z=Vertex[6].z=(float)(pos.z+orient.c0.z+orient.c1.z-orient.c2.z);
	
	Vertex[23].x=Vertex[14].x=Vertex[7].x=(float)(pos.x-orient.c0.x+orient.c1.x-orient.c2.x);
	Vertex[23].y=Vertex[14].y=Vertex[7].y=(float)(pos.y-orient.c0.y+orient.c1.y-orient.c2.y);
	Vertex[23].z=Vertex[14].z=Vertex[7].z=(float)(pos.z-orient.c0.z+orient.c1.z-orient.c2.z);
}

WORD Idx[36]=
{ 0,  3,  1,
  1,  3,  2,
  4,  5,  7,
  5,  6,  7,
  9, 10,  8,
  8, 10, 11,
 12, 15, 13,
 14, 13, 15,
 18, 16, 17,
 16, 18, 19,
 23, 21, 20,
 21, 23, 22
};

MyCube **c=NULL;

VOID AppOutputText( LPDIRECT3DDEVICE3 pd3dDevice, DWORD x, DWORD y, CHAR* str )
{
    LPDIRECTDRAWSURFACE4 pddsRenderSurface;
    if( FAILED( pd3dDevice->GetRenderTarget( &pddsRenderSurface ) ) )
        return;

    // Get a DC for the surface. Then, write out the buffer
    HDC hDC;
    if( SUCCEEDED( pddsRenderSurface->GetDC(&hDC) ) )
    {
        SetTextColor( hDC, RGB(255,255,0) );
        SetBkMode( hDC, TRANSPARENT );
        ExtTextOut( hDC, x, y, 0, NULL, str, strlen(str), NULL );
    
        pddsRenderSurface->ReleaseDC(hDC);
    }
    pddsRenderSurface->Release();
}

VOID AppShowStats()
{
    static FLOAT fFPS      = 0.0f;
    static FLOAT fLastTime = 0.0f;
    static DWORD dwFrames  = 0L;

	// Keep track of the time lapse and frame count
	FLOAT fTime = timeGetTime() * 0.001f; // Get current time in seconds
	++dwFrames;

	// Update the frame rate once per second
	if( fTime - fLastTime > 1.0f )
    {
        fFPS      = dwFrames / (fTime - fLastTime);
        fLastTime = fTime;
        dwFrames  = 0L;
    }

    // Setup the text buffer to write out
    CHAR buffer[80];
    sprintf( buffer, "%7.02f fps", fFPS );
    AppOutputText( g_pd3dDevice, 0, 0, buffer );
}

HRESULT App_Render( LPDIRECT3DDEVICE3 pd3dDevice, 
				   LPDIRECT3DVIEWPORT3 pvViewport,
				   D3DRECT* prcViewportRect )
{
	pvViewport->Clear2( 1UL, prcViewportRect, D3DCLEAR_TARGET|D3DCLEAR_ZBUFFER,
		0x00000000,
		1.0f, 0L );
	
	if( FAILED( pd3dDevice->BeginScene() ) )
		return E_FAIL;
	
	if (c) {
		int i;
		for (i=0;i<NCUBES;i++){
			SetVertices(c[i]->pos, c[i]->orient);
			pd3dDevice->DrawIndexedPrimitive(D3DPT_TRIANGLELIST,VertexType,Vertex,24,Idx,36,0);
		}
	}
	
	
	pd3dDevice->EndScene();

    AppShowStats();

	return S_OK;
}

boolean HandleEvents(){
	BOOL bGotMsg,changed=TRUE;
	PeekMessage( &msg, NULL, 0U, 0U, PM_NOREMOVE );
	g_bReady = TRUE;
	
    while( WM_QUIT != msg.message  )
    {
		if( g_bActive )
			bGotMsg = PeekMessage( &msg, NULL, 0U, 0U, PM_REMOVE );
		else
			bGotMsg = GetMessage( &msg, NULL, 0U, 0U );
		
		if( bGotMsg )
        {
            if( 0 == TranslateAccelerator( MainWin, hAccel, &msg ) )
			{
				TranslateMessage( &msg );
				DispatchMessage( &msg );
			}
			changed=TRUE;
        }
		else
		{
			if( g_bActive && g_bReady )
				if (changed) {
					Render3DEnvironment();
					changed=FALSE;
					return true;
				}
		}
    }
	return false;
}

boolean RenderCubes(MyCube* cube[]){
	c=cube;
	return HandleEvents();
}

