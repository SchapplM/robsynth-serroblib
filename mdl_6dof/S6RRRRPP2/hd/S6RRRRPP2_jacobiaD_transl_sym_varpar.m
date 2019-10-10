% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(7) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:16
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
	t38 = qJD(2) + qJD(3);
	t39 = qJ(2) + qJ(3);
	t37 = cos(t39);
	t59 = r_i_i_C(2) * t37;
	t36 = sin(t39);
	t61 = r_i_i_C(1) * t36;
	t49 = t59 + t61;
	t47 = t49 * t38;
	t40 = sin(qJ(2));
	t62 = pkin(2) * t40;
	t63 = qJD(2) * t62 + t47;
	t60 = r_i_i_C(2) * t36;
	t58 = r_i_i_C(3) + pkin(8) + pkin(7);
	t57 = t37 * t38;
	t41 = sin(qJ(1));
	t56 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t55 = qJD(1) * t43;
	t42 = cos(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(1) * t57;
	t52 = t38 * t60;
	t51 = qJD(1) * t59;
	t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
	t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
	t31 = t41 * t52;
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0, 0; 0, -t63, -t47, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:17
	% EndTime: 2019-10-10 12:22:17
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (270->58), mult. (394->95), div. (0->0), fcn. (297->8), ass. (0->54)
	t258 = qJ(2) + qJ(3);
	t255 = sin(t258);
	t262 = cos(qJ(4));
	t310 = r_i_i_C(1) * t262 + pkin(3);
	t313 = t255 * t310;
	t259 = sin(qJ(4));
	t293 = qJD(4) * t262;
	t256 = cos(t258);
	t257 = qJD(2) + qJD(3);
	t301 = t256 * t257;
	t312 = t255 * t293 + t259 * t301;
	t307 = pkin(9) + r_i_i_C(3);
	t288 = t307 * t256;
	t260 = sin(qJ(2));
	t302 = pkin(2) * qJD(2);
	t290 = t260 * t302;
	t305 = pkin(3) * t255;
	t311 = -t290 + (t288 - t305) * t257;
	t294 = qJD(4) * t259;
	t281 = t255 * t294;
	t308 = r_i_i_C(1) * t281 + t312 * r_i_i_C(2);
	t306 = pkin(2) * t260;
	t303 = r_i_i_C(2) * t259;
	t261 = sin(qJ(1));
	t300 = t257 * t261;
	t299 = t257 * t262;
	t264 = cos(qJ(1));
	t298 = t257 * t264;
	t297 = t262 * t264;
	t296 = qJD(1) * t261;
	t295 = qJD(1) * t264;
	t292 = t255 * t303;
	t291 = qJD(1) * t303;
	t289 = t307 * t255;
	t287 = t307 * t261;
	t286 = t255 * t299;
	t276 = qJD(4) * t256 - qJD(1);
	t275 = qJD(1) * t256 - qJD(4);
	t274 = t310 * t256;
	t273 = t310 * t264;
	t272 = t308 * t264 + t296 * t313;
	t271 = t276 * t259;
	t270 = t264 * t255 * t291 + t308 * t261 + t295 * t288;
	t263 = cos(qJ(2));
	t269 = -t263 * pkin(2) - pkin(3) * t256 - pkin(1) - t289;
	t268 = t255 * t298 + t275 * t261;
	t267 = -t263 * t302 + (-t274 - t289) * t257;
	t266 = -t256 * r_i_i_C(2) * t293 + (-t256 * t294 - t286) * r_i_i_C(1) + t307 * t301 + (-t305 + t292) * t257;
	t265 = -pkin(8) - pkin(7);
	t239 = -t275 * t297 + (t271 + t286) * t261;
	t238 = t276 * t262 * t261 + (-t255 * t300 + t275 * t264) * t259;
	t237 = t268 * t262 + t264 * t271;
	t236 = t268 * t259 - t276 * t297;
	t1 = [t239 * r_i_i_C(1) + t238 * r_i_i_C(2) - t311 * t261 + (t261 * t265 + t269 * t264) * qJD(1), (-t288 - t292 + t306) * t296 + t267 * t264 + t272, (-t261 * t291 - t307 * t298) * t255 + (-qJD(1) * t287 - t257 * t273) * t256 + t272, t236 * r_i_i_C(1) + t237 * r_i_i_C(2), 0, 0; -t237 * r_i_i_C(1) + t236 * r_i_i_C(2) + t311 * t264 + (t269 * t261 - t264 * t265) * qJD(1), (-t306 - t313) * t295 + t267 * t261 + t270, -t274 * t300 + (-qJD(1) * t273 - t257 * t287) * t255 + t270, -t238 * r_i_i_C(1) + t239 * r_i_i_C(2), 0, 0; 0, t266 - t290, t266, (-t256 * t299 + t281) * r_i_i_C(2) - t312 * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:17
	% EndTime: 2019-10-10 12:22:18
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (474->71), mult. (702->99), div. (0->0), fcn. (565->8), ass. (0->56)
	t316 = sin(qJ(4));
	t319 = cos(qJ(4));
	t321 = cos(qJ(1));
	t357 = qJD(4) * t321;
	t318 = sin(qJ(1));
	t360 = qJD(1) * t318;
	t330 = t316 * t357 + t319 * t360;
	t369 = r_i_i_C(1) + pkin(4);
	t367 = r_i_i_C(3) + qJ(5);
	t374 = t367 * t316;
	t315 = qJ(2) + qJ(3);
	t313 = cos(t315);
	t370 = pkin(9) + r_i_i_C(2);
	t351 = t370 * t313;
	t314 = qJD(2) + qJD(3);
	t350 = t370 * t314;
	t317 = sin(qJ(2));
	t366 = pkin(2) * qJD(2);
	t354 = t317 * t366;
	t356 = qJD(5) * t316;
	t312 = sin(t315);
	t365 = t312 * t314;
	t373 = -pkin(3) * t365 + (t350 + t356) * t313 - t354;
	t358 = qJD(4) * t319;
	t328 = -t367 * t358 - t356;
	t368 = pkin(2) * t317;
	t364 = t313 * t314;
	t363 = t313 * t318;
	t362 = t316 * t318;
	t361 = t319 * t321;
	t359 = qJD(1) * t321;
	t355 = t319 * qJD(5);
	t353 = t369 * t316;
	t352 = t370 * t312;
	t349 = t318 * t365;
	t348 = t321 * t365;
	t343 = qJD(4) * t362;
	t341 = t319 * t357;
	t340 = t369 * t312 * t343 + t359 * t351;
	t335 = t313 * t361 + t362;
	t320 = cos(qJ(2));
	t333 = -pkin(2) * t320 - pkin(3) * t313 - pkin(1) - t352;
	t332 = (t369 * t330 + (pkin(3) + t374) * t360) * t312;
	t331 = -t369 * t319 - t374;
	t329 = t316 * t359 + t318 * t358;
	t327 = -pkin(3) + t331;
	t326 = t327 * t314;
	t325 = t312 * t326 + t370 * t364 + (-qJD(4) * t353 - t328) * t313;
	t324 = t328 * t312 + (t327 * t313 - t352) * t314;
	t323 = -t320 * t366 + t324;
	t322 = -pkin(8) - pkin(7);
	t284 = t335 * qJD(1) - t313 * t343 - t319 * t349 - t341;
	t283 = t329 * t313 - t316 * t349 - t330;
	t282 = t330 * t313 + t319 * t348 - t329;
	t281 = t316 * t348 - t313 * t341 - t343 + (t313 * t362 + t361) * qJD(1);
	t1 = [-t321 * t355 - t369 * t284 - t367 * t283 - t373 * t318 + (t318 * t322 + t333 * t321) * qJD(1), (-t351 + t368) * t360 + t323 * t321 + t332, t324 * t321 - t351 * t360 + t332, t335 * qJD(5) + t369 * t281 - t367 * t282, -t281, 0; -t318 * t355 - t369 * t282 - t367 * t281 + t373 * t321 + (t333 * t318 - t321 * t322) * qJD(1), (t327 * t312 - t368) * t359 + t323 * t318 + t340, t326 * t363 + ((-t350 + t328) * t318 + t327 * t359) * t312 + t340, -(t321 * t316 - t319 * t363) * qJD(5) + t367 * t284 - t369 * t283, t283, 0; 0, t325 - t354, t325, (t367 * t319 - t353) * t364 + (t331 * qJD(4) + t355) * t312, t312 * t358 + t316 * t364, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:17
	% EndTime: 2019-10-10 12:22:17
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (616->78), mult. (887->104), div. (0->0), fcn. (714->8), ass. (0->57)
	t334 = r_i_i_C(3) + qJ(6);
	t275 = sin(qJ(4));
	t326 = r_i_i_C(2) + qJ(5);
	t336 = t326 * t275;
	t278 = cos(qJ(4));
	t280 = cos(qJ(1));
	t315 = qJD(4) * t280;
	t277 = sin(qJ(1));
	t318 = qJD(1) * t277;
	t291 = t275 * t315 + t278 * t318;
	t314 = qJD(5) * t275;
	t316 = qJD(4) * t278;
	t335 = -t326 * t316 - t314;
	t312 = r_i_i_C(1) + pkin(5) + pkin(4);
	t274 = qJ(2) + qJ(3);
	t271 = sin(t274);
	t272 = cos(t274);
	t273 = qJD(2) + qJD(3);
	t276 = sin(qJ(2));
	t324 = pkin(2) * qJD(2);
	t310 = t276 * t324;
	t311 = pkin(9) - t334;
	t332 = (t311 * t273 + t314) * t272 - (pkin(3) * t273 + qJD(6)) * t271 - qJD(1) * (-pkin(8) - pkin(7)) - t310;
	t331 = t271 * t312;
	t328 = pkin(2) * t276;
	t327 = pkin(9) * t272;
	t323 = t272 * t273;
	t322 = t273 * t277;
	t321 = t273 * t280;
	t320 = t277 * t275;
	t319 = t280 * t278;
	t317 = qJD(1) * t280;
	t313 = t278 * qJD(5);
	t309 = t271 * t322;
	t308 = t271 * t321;
	t307 = t271 * t318;
	t306 = t272 * t318;
	t303 = qJD(4) * t320;
	t301 = t278 * t315;
	t299 = t312 * t275;
	t294 = t303 * t331 + t334 * t309 + t317 * t327;
	t292 = t272 * t319 + t320;
	t290 = t275 * t317 + t277 * t316;
	t289 = -t312 * t278 - t336;
	t288 = -pkin(3) + t289;
	t279 = cos(qJ(2));
	t287 = -t313 + (-t279 * pkin(2) - pkin(3) * t272 - t311 * t271 - pkin(1)) * qJD(1);
	t286 = t334 * (t306 + t308) + t291 * t331 + (pkin(3) + t336) * t307;
	t285 = t288 * t271 - t272 * t334;
	t284 = (-pkin(9) * t273 + t335) * t271 + (t288 * t273 - qJD(6)) * t272;
	t283 = -t279 * t324 + t284;
	t282 = pkin(9) * t323 - t271 * qJD(6) + t285 * t273 + (-qJD(4) * t299 - t335) * t272;
	t236 = t292 * qJD(1) - t272 * t303 - t278 * t309 - t301;
	t235 = t290 * t272 - t275 * t309 - t291;
	t234 = t291 * t272 + t278 * t308 - t290;
	t233 = t275 * t308 - t272 * t301 - t303 + (t272 * t320 + t319) * qJD(1);
	t1 = [-t326 * t235 - t312 * t236 - t332 * t277 + t287 * t280, t286 + t283 * t280 + (-t327 + t328) * t318, -pkin(9) * t306 + t284 * t280 + t286, t292 * qJD(5) + t312 * t233 - t326 * t234, -t233, -t272 * t321 + t307; -t326 * t233 - t312 * t234 + t287 * t277 + t332 * t280, (t285 - t328) * t317 + t283 * t277 + t294, t284 * t277 + t285 * t317 + t294, -(-t277 * t272 * t278 + t280 * t275) * qJD(5) + t326 * t236 - t312 * t235, t235, -t271 * t317 - t272 * t322; 0, t282 - t310, t282, (t326 * t278 - t299) * t323 + (t289 * qJD(4) + t313) * t271, t271 * t316 + t275 * t323, -t273 * t271;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end