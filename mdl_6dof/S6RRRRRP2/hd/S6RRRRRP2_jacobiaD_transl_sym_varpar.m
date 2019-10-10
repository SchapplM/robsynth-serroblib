% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRP2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:58
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
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
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
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
	% StartTime: 2019-10-10 12:58:34
	% EndTime: 2019-10-10 12:58:34
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
	% StartTime: 2019-10-10 12:58:35
	% EndTime: 2019-10-10 12:58:36
	% DurationCPUTime: 0.54s
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
	t287 = t278 * t264 * t309 + t329 * t275 + t306 * t315;
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
	% StartTime: 2019-10-10 12:58:36
	% EndTime: 2019-10-10 12:58:36
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (953->84), mult. (900->108), div. (0->0), fcn. (706->10), ass. (0->66)
	t330 = sin(qJ(5));
	t333 = cos(qJ(5));
	t335 = cos(qJ(1));
	t375 = qJD(5) * t335;
	t332 = sin(qJ(1));
	t378 = qJD(1) * t332;
	t346 = t330 * t375 + t333 * t378;
	t388 = -r_i_i_C(1) - pkin(5);
	t385 = r_i_i_C(3) + qJ(6);
	t393 = t385 * t330;
	t329 = qJ(2) + qJ(3);
	t326 = qJ(4) + t329;
	t322 = cos(t326);
	t389 = pkin(10) + r_i_i_C(2);
	t368 = t389 * t322;
	t327 = qJD(2) + qJD(3);
	t323 = qJD(4) + t327;
	t367 = t389 * t323;
	t331 = sin(qJ(2));
	t384 = pkin(2) * qJD(2);
	t370 = t331 * t384;
	t324 = sin(t329);
	t386 = pkin(3) * t327;
	t372 = t324 * t386;
	t374 = qJD(6) * t330;
	t321 = sin(t326);
	t383 = t321 * t323;
	t392 = -pkin(4) * t383 + (t367 + t374) * t322 - t370 - t372;
	t376 = qJD(5) * t333;
	t344 = -t385 * t376 - t374;
	t325 = cos(t329);
	t347 = t388 * t333 - t393;
	t343 = -pkin(4) + t347;
	t369 = t389 * t321;
	t338 = t344 * t321 + (t343 * t322 - t369) * t323;
	t336 = -t325 * t386 + t338;
	t387 = pkin(3) * t324;
	t382 = t322 * t323;
	t381 = t322 * t332;
	t380 = t332 * t330;
	t379 = t335 * t333;
	t377 = qJD(1) * t335;
	t373 = t333 * qJD(6);
	t366 = t388 * t330;
	t365 = t332 * t383;
	t364 = t335 * t383;
	t359 = qJD(5) * t380;
	t357 = t333 * t375;
	t356 = -t388 * t321 * t359 + t377 * t368;
	t350 = t322 * t379 + t380;
	t334 = cos(qJ(2));
	t349 = -t334 * pkin(2) - pkin(3) * t325 - pkin(4) * t322 - pkin(1) - t369;
	t348 = (-t388 * t346 + (pkin(4) + t393) * t378) * t321;
	t345 = t330 * t377 + t332 * t376;
	t342 = t343 * t321;
	t341 = t343 * t323;
	t340 = t321 * t341 + t389 * t382 + (qJD(5) * t366 - t344) * t322;
	t339 = t340 - t372;
	t337 = -t334 * t384 + t336;
	t328 = -pkin(9) - pkin(8) - pkin(7);
	t313 = -t331 * pkin(2) - t387;
	t290 = t350 * qJD(1) - t322 * t359 - t333 * t365 - t357;
	t289 = t345 * t322 - t330 * t365 - t346;
	t288 = t346 * t322 + t333 * t364 - t345;
	t287 = t330 * t364 - t322 * t357 - t359 + (t322 * t380 + t379) * qJD(1);
	t1 = [-t335 * t373 + t388 * t290 - t385 * t289 - t392 * t332 + (t332 * t328 + t349 * t335) * qJD(1), (-t313 - t368) * t378 + t337 * t335 + t348, (-t368 + t387) * t378 + t336 * t335 + t348, t338 * t335 - t368 * t378 + t348, t350 * qJD(6) - t388 * t287 - t385 * t288, -t287; -t332 * t373 + t388 * t288 - t385 * t287 + t392 * t335 + (-t335 * t328 + t349 * t332) * qJD(1), (t313 + t342) * t377 + t337 * t332 + t356, (t342 - t387) * t377 + t336 * t332 + t356, t341 * t381 + ((-t367 + t344) * t332 + t343 * t377) * t321 + t356, -(t335 * t330 - t333 * t381) * qJD(6) + t385 * t290 + t388 * t289, t289; 0, t339 - t370, t339, t340, (t385 * t333 + t366) * t382 + (t347 * qJD(5) + t373) * t321, t321 * t376 + t330 * t382;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end