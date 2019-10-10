% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
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
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:03
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:03
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-10 13:18:03
	% EndTime: 2019-10-10 13:18:03
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (202->36), mult. (170->47), div. (0->0), fcn. (109->8), ass. (0->37)
	t54 = sin(qJ(2));
	t51 = qJD(2) + qJD(3);
	t47 = qJD(4) + t51;
	t53 = qJ(2) + qJ(3);
	t50 = qJ(4) + t53;
	t46 = cos(t50);
	t77 = r_i_i_C(2) * t46;
	t45 = sin(t50);
	t79 = r_i_i_C(1) * t45;
	t62 = t77 + t79;
	t60 = t62 * t47;
	t48 = sin(t53);
	t80 = pkin(3) * t48;
	t58 = -t51 * t80 - t60;
	t73 = pkin(2) * qJD(2);
	t82 = -t54 * t73 + t58;
	t75 = t46 * t47;
	t68 = r_i_i_C(1) * t75;
	t49 = cos(t53);
	t74 = t49 * t51;
	t81 = -pkin(3) * t74 - t68;
	t78 = r_i_i_C(2) * t45;
	t76 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7);
	t55 = sin(qJ(1));
	t72 = qJD(1) * t55;
	t57 = cos(qJ(1));
	t71 = qJD(1) * t57;
	t67 = t47 * t78;
	t65 = qJD(1) * t77;
	t66 = t55 * t65 + t57 * t67 + t72 * t79;
	t56 = cos(qJ(2));
	t63 = -t56 * t73 + t81;
	t61 = -pkin(2) * t56 - pkin(3) * t49 - r_i_i_C(1) * t46 - pkin(1) + t78;
	t59 = -t57 * t68 + t66;
	t44 = -pkin(2) * t54 - t80;
	t39 = t55 * t67;
	t1 = [-t82 * t55 + (-t76 * t55 + t61 * t57) * qJD(1), -t44 * t72 + t63 * t57 + t66, (t48 * t72 - t57 * t74) * pkin(3) + t59, t59, 0, 0; t82 * t57 + (t61 * t55 + t76 * t57) * qJD(1), t39 + t63 * t55 + (t44 - t62) * t71, t39 + t81 * t55 + (-t62 - t80) * t71, -t57 * t65 + t39 + (-t45 * t71 - t55 * t75) * r_i_i_C(1), 0, 0; 0, t82, t58, -t60, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:04
	% EndTime: 2019-10-10 13:18:04
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (576->71), mult. (524->104), div. (0->0), fcn. (386->10), ass. (0->64)
	t272 = qJ(2) + qJ(3);
	t269 = qJ(4) + t272;
	t264 = sin(t269);
	t276 = cos(qJ(5));
	t331 = r_i_i_C(1) * t276 + pkin(4);
	t292 = t331 * t264;
	t273 = sin(qJ(5));
	t313 = qJD(5) * t276;
	t265 = cos(t269);
	t270 = qJD(2) + qJD(3);
	t266 = qJD(4) + t270;
	t321 = t265 * t266;
	t333 = t264 * t313 + t273 * t321;
	t328 = pkin(10) + r_i_i_C(3);
	t306 = t328 * t265;
	t274 = sin(qJ(2));
	t322 = pkin(2) * qJD(2);
	t308 = t274 * t322;
	t267 = sin(t272);
	t326 = pkin(3) * t270;
	t312 = t267 * t326;
	t325 = pkin(4) * t264;
	t332 = (t306 - t325) * t266 - t308 - t312;
	t314 = qJD(5) * t273;
	t299 = t264 * t314;
	t329 = r_i_i_C(1) * t299 + t333 * r_i_i_C(2);
	t268 = cos(t272);
	t291 = t331 * t265;
	t307 = t328 * t264;
	t281 = (-t291 - t307) * t266 - t268 * t326;
	t327 = pkin(3) * t267;
	t323 = r_i_i_C(2) * t273;
	t275 = sin(qJ(1));
	t320 = t266 * t275;
	t319 = t266 * t276;
	t278 = cos(qJ(1));
	t318 = t266 * t278;
	t317 = t276 * t278;
	t316 = qJD(1) * t275;
	t315 = qJD(1) * t278;
	t310 = t264 * t323;
	t309 = qJD(1) * t323;
	t305 = t328 * t275;
	t304 = t264 * t319;
	t294 = qJD(5) * t265 - qJD(1);
	t293 = qJD(1) * t265 - qJD(5);
	t290 = t331 * t278;
	t289 = t329 * t278 + t316 * t292;
	t288 = t294 * t273;
	t287 = t278 * t264 * t309 + t329 * t275 + t315 * t306;
	t286 = -t306 - t310;
	t277 = cos(qJ(2));
	t285 = -t277 * pkin(2) - pkin(3) * t268 - pkin(4) * t265 - pkin(1) - t307;
	t284 = t264 * t318 + t293 * t275;
	t282 = -t277 * t322 + t281;
	t280 = -t265 * r_i_i_C(2) * t313 + (-t265 * t314 - t304) * r_i_i_C(1) + t328 * t321 + (-t325 + t310) * t266;
	t279 = t280 - t312;
	t271 = -pkin(9) - pkin(8) - pkin(7);
	t263 = -t274 * pkin(2) - t327;
	t245 = -t293 * t317 + (t288 + t304) * t275;
	t244 = t294 * t276 * t275 + (-t264 * t320 + t293 * t278) * t273;
	t243 = t284 * t276 + t278 * t288;
	t242 = t284 * t273 - t294 * t317;
	t1 = [t245 * r_i_i_C(1) + t244 * r_i_i_C(2) - t332 * t275 + (t271 * t275 + t285 * t278) * qJD(1), (-t263 + t286) * t316 + t282 * t278 + t289, (t286 + t327) * t316 + t281 * t278 + t289, (-t275 * t309 - t328 * t318) * t264 + (-qJD(1) * t305 - t266 * t290) * t265 + t289, t242 * r_i_i_C(1) + t243 * r_i_i_C(2), 0; -t243 * r_i_i_C(1) + t242 * r_i_i_C(2) + t332 * t278 + (-t271 * t278 + t285 * t275) * qJD(1), (t263 - t292) * t315 + t282 * t275 + t287, (-t292 - t327) * t315 + t281 * t275 + t287, -t291 * t320 + (-qJD(1) * t290 - t266 * t305) * t264 + t287, -t244 * r_i_i_C(1) + t245 * r_i_i_C(2), 0; 0, t279 - t308, t279, t280, (-t265 * t319 + t299) * r_i_i_C(2) - t333 * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:04
	% EndTime: 2019-10-10 13:18:05
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (912->91), mult. (706->119), div. (0->0), fcn. (533->12), ass. (0->79)
	t307 = qJ(2) + qJ(3);
	t302 = qJ(4) + t307;
	t294 = sin(t302);
	t306 = qJ(5) + qJ(6);
	t298 = sin(t306);
	t300 = cos(t306);
	t303 = qJD(5) + qJD(6);
	t359 = t300 * t303;
	t295 = cos(t302);
	t304 = qJD(2) + qJD(3);
	t297 = qJD(4) + t304;
	t363 = t295 * t297;
	t384 = t294 * t359 + t298 * t363;
	t314 = -pkin(11) - pkin(10);
	t308 = sin(qJ(5));
	t366 = pkin(5) * qJD(5);
	t350 = t308 * t366;
	t383 = t297 * t314 + t350;
	t311 = cos(qJ(5));
	t296 = t311 * pkin(5) + pkin(4);
	t381 = r_i_i_C(1) * t300 + t296;
	t382 = t294 * t381 + t295 * t314;
	t309 = sin(qJ(2));
	t367 = pkin(2) * qJD(2);
	t347 = t309 * t367;
	t299 = sin(t307);
	t373 = pkin(3) * t304;
	t353 = t299 * t373;
	t365 = t294 * t297;
	t368 = r_i_i_C(3) - t314;
	t372 = pkin(5) * t308;
	t380 = (t368 * t297 - t350) * t295 + (pkin(9) + pkin(8) + pkin(7) + t372) * qJD(1) - t296 * t365 - t347 - t353;
	t364 = t294 * t298;
	t343 = t303 * t364;
	t379 = r_i_i_C(1) * t343 + t384 * r_i_i_C(2) + t383 * t294;
	t313 = cos(qJ(1));
	t334 = t295 * t303 - qJD(1);
	t378 = t313 * t334;
	t301 = cos(t307);
	t321 = (-r_i_i_C(3) * t294 - t295 * t381) * t297;
	t318 = -t301 * t373 + t321;
	t356 = qJD(1) * t295;
	t333 = -t303 + t356;
	t310 = sin(qJ(1));
	t346 = t310 * t365;
	t376 = t333 * t313 - t346;
	t374 = pkin(3) * t299;
	t370 = r_i_i_C(2) * t300;
	t369 = r_i_i_C(3) * t295;
	t361 = t297 * t313;
	t345 = t294 * t361;
	t320 = t333 * t310 + t345;
	t268 = t320 * t298 - t300 * t378;
	t269 = t298 * t378 + t320 * t300;
	t358 = t268 * r_i_i_C(1) + t269 * r_i_i_C(2);
	t327 = t334 * t310;
	t270 = t376 * t298 + t300 * t327;
	t271 = t298 * t327 - t376 * t300;
	t357 = -t270 * r_i_i_C(1) + t271 * r_i_i_C(2);
	t355 = qJD(1) * t310;
	t354 = qJD(1) * t313;
	t351 = r_i_i_C(2) * t364;
	t349 = t311 * t366;
	t348 = r_i_i_C(2) * qJD(1) * t298;
	t331 = -qJD(5) + t356;
	t330 = -r_i_i_C(1) * t298 - t370;
	t328 = t381 * t297;
	t326 = (-qJD(5) * t295 + qJD(1)) * t311;
	t325 = -t351 - t369;
	t324 = t313 * t294 * t348 + t310 * t379 + t354 * t369;
	t322 = t379 * t313 + t382 * t355;
	t312 = cos(qJ(2));
	t319 = -t312 * t367 + t318;
	t317 = t349 + (-t312 * pkin(2) - pkin(3) * t301 - t368 * t294 - t295 * t296 - pkin(1)) * qJD(1);
	t316 = t297 * t351 + r_i_i_C(3) * t363 - t294 * t328 + (t330 * t303 - t383) * t295;
	t315 = t316 - t353;
	t293 = -t309 * pkin(2) - t374;
	t286 = r_i_i_C(2) * t343;
	t1 = [t271 * r_i_i_C(1) + t270 * r_i_i_C(2) - t310 * t380 + t317 * t313, (-t293 + t325) * t355 + t319 * t313 + t322, (t325 + t374) * t355 + t318 * t313 + t322, (-r_i_i_C(3) * t361 - t310 * t348) * t294 + (-r_i_i_C(3) * t355 - t313 * t328) * t295 + t322, (t313 * t326 + (t331 * t310 + t345) * t308) * pkin(5) + t358, t358; -t269 * r_i_i_C(1) + t268 * r_i_i_C(2) + t317 * t310 + t313 * t380, t319 * t310 + (t293 - t382) * t354 + t324, t318 * t310 + (-t382 - t374) * t354 + t324, t310 * t321 - t354 * t382 + t324, (t310 * t326 + (-t331 * t313 + t346) * t308) * pkin(5) + t357, t357; 0, t315 - t347, t315, t316, t286 + (-r_i_i_C(1) * t359 - t349) * t294 + (t330 - t372) * t363, -t384 * r_i_i_C(1) - t363 * t370 + t286;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end