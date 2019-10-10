% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:45
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
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
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(10);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (43->18), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t24 = sin(qJ(3));
	t25 = cos(qJ(3));
	t26 = (r_i_i_C(1) * t24 + r_i_i_C(2) * t25) * qJD(3);
	t33 = pkin(7) + r_i_i_C(3);
	t32 = qJD(1) * t24;
	t31 = qJD(1) * t25;
	t30 = qJD(3) * t24;
	t29 = qJD(3) * t25;
	t27 = -r_i_i_C(1) * t25 + r_i_i_C(2) * t24 - pkin(2);
	t23 = qJ(1) + pkin(10);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (119->29), mult. (118->40), div. (0->0), fcn. (75->8), ass. (0->26)
	t44 = qJD(3) + qJD(4);
	t46 = qJ(3) + qJ(4);
	t42 = sin(t46);
	t43 = cos(t46);
	t63 = r_i_i_C(2) * t43;
	t54 = r_i_i_C(1) * t42 + t63;
	t52 = t54 * t44;
	t47 = sin(qJ(3));
	t65 = pkin(3) * t47;
	t66 = qJD(3) * t65 + t52;
	t64 = r_i_i_C(2) * t42;
	t62 = r_i_i_C(3) + pkin(8) + pkin(7);
	t61 = t43 * t44;
	t60 = qJD(1) * t42;
	t48 = cos(qJ(3));
	t59 = qJD(3) * t48;
	t58 = r_i_i_C(1) * t61;
	t57 = t44 * t64;
	t56 = qJD(1) * t63;
	t53 = -t48 * pkin(3) - r_i_i_C(1) * t43 - pkin(2) + t64;
	t45 = qJ(1) + pkin(10);
	t40 = sin(t45);
	t41 = cos(t45);
	t51 = (t57 - t58) * t41 + (t60 * r_i_i_C(1) + t56) * t40;
	t35 = t40 * t57;
	t1 = [t66 * t40 + (-cos(qJ(1)) * pkin(1) - t62 * t40 + t53 * t41) * qJD(1), 0, (qJD(1) * t40 * t47 - t41 * t59) * pkin(3) + t51, t51, 0, 0; -t66 * t41 + (-sin(qJ(1)) * pkin(1) + t62 * t41 + t53 * t40) * qJD(1), 0, t35 + (-pkin(3) * t59 - t58) * t40 + (-t54 - t65) * t41 * qJD(1), -t41 * t56 + t35 + (-t40 * t61 - t41 * t60) * r_i_i_C(1), 0, 0; 0, 0, -t66, -t52, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:25
	% EndTime: 2019-10-10 01:45:25
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (382->56), mult. (398->81), div. (0->0), fcn. (299->10), ass. (0->53)
	t265 = qJ(3) + qJ(4);
	t261 = sin(t265);
	t268 = cos(qJ(5));
	t315 = r_i_i_C(1) * t268 + pkin(4);
	t280 = t315 * t261;
	t262 = cos(t265);
	t298 = qJD(5) * t268;
	t263 = qJD(3) + qJD(4);
	t266 = sin(qJ(5));
	t304 = t263 * t266;
	t318 = t261 * t298 + t262 * t304;
	t310 = pkin(9) + r_i_i_C(3);
	t294 = t310 * t262;
	t317 = (-pkin(4) * t261 + t294) * t263;
	t267 = sin(qJ(3));
	t306 = pkin(3) * qJD(3);
	t296 = t267 * t306;
	t316 = -t296 + t317;
	t299 = qJD(5) * t266;
	t287 = t261 * t299;
	t313 = r_i_i_C(1) * t287 + t318 * r_i_i_C(2);
	t300 = qJD(1) * t262;
	t281 = -qJD(5) + t300;
	t312 = t268 * t281;
	t282 = qJD(5) * t262 - qJD(1);
	t293 = t261 * t304;
	t311 = t282 * t268 - t293;
	t309 = pkin(3) * t267;
	t303 = t263 * t268;
	t264 = qJ(1) + pkin(10);
	t259 = sin(t264);
	t302 = qJD(1) * t259;
	t260 = cos(t264);
	t301 = qJD(1) * t260;
	t297 = r_i_i_C(2) * t261 * t266;
	t295 = t310 * t261;
	t292 = t261 * t303;
	t279 = t313 * t260 + t302 * t280;
	t278 = t281 * t266;
	t277 = t310 * t260 * t300 + t313 * t259 + t297 * t301;
	t276 = -t294 - t297;
	t269 = cos(qJ(3));
	t275 = -t269 * pkin(3) - pkin(4) * t262 - pkin(2) - t295;
	t274 = t282 * t266 + t292;
	t273 = (-t262 * t315 - t295) * t263;
	t272 = -t269 * t306 + t273;
	t271 = (-t262 * t299 - t292) * r_i_i_C(1) + (-t262 * t298 + t293) * r_i_i_C(2) + t317;
	t270 = -pkin(8) - pkin(7);
	t243 = t274 * t259 - t260 * t312;
	t242 = t311 * t259 + t260 * t278;
	t241 = t259 * t312 + t274 * t260;
	t240 = t259 * t278 - t311 * t260;
	t1 = [t243 * r_i_i_C(1) + t242 * r_i_i_C(2) - t316 * t259 + (-cos(qJ(1)) * pkin(1) + t259 * t270 + t275 * t260) * qJD(1), 0, (t276 + t309) * t302 + t272 * t260 + t279, t260 * t273 + t276 * t302 + t279, t240 * r_i_i_C(1) + t241 * r_i_i_C(2), 0; -t241 * r_i_i_C(1) + t240 * r_i_i_C(2) + t316 * t260 + (-sin(qJ(1)) * pkin(1) - t260 * t270 + t275 * t259) * qJD(1), 0, (-t280 - t309) * t301 + t272 * t259 + t277, t259 * t273 - t280 * t301 + t277, -t242 * r_i_i_C(1) + t243 * r_i_i_C(2), 0; 0, 0, t271 - t296, t271, (-t262 * t303 + t287) * r_i_i_C(2) - t318 * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:25
	% EndTime: 2019-10-10 01:45:26
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (672->73), mult. (706->100), div. (0->0), fcn. (567->10), ass. (0->57)
	t313 = qJ(1) + pkin(10);
	t309 = cos(t313);
	t317 = cos(qJ(5));
	t308 = sin(t313);
	t354 = qJD(1) * t308;
	t342 = t317 * t354;
	t315 = sin(qJ(5));
	t352 = qJD(5) * t315;
	t370 = t309 * t352 + t342;
	t363 = r_i_i_C(1) + pkin(5);
	t361 = r_i_i_C(3) + qJ(6);
	t369 = t361 * t315;
	t314 = qJ(3) + qJ(4);
	t311 = cos(t314);
	t364 = pkin(9) + r_i_i_C(2);
	t346 = t364 * t311;
	t350 = qJD(6) * t315;
	t351 = qJD(5) * t317;
	t368 = t361 * t351 + t350;
	t312 = qJD(3) + qJD(4);
	t316 = sin(qJ(3));
	t360 = pkin(3) * qJD(3);
	t348 = t316 * t360;
	t310 = sin(t314);
	t358 = t310 * t312;
	t366 = -pkin(4) * t358 + (t364 * t312 + t350) * t311 - t348;
	t362 = pkin(3) * t316;
	t359 = t308 * t315;
	t357 = t311 * t312;
	t356 = t311 * t317;
	t355 = t312 * t315;
	t353 = qJD(1) * t309;
	t349 = t317 * qJD(6);
	t347 = t364 * t310;
	t345 = t310 * t355;
	t344 = t317 * t358;
	t340 = t311 * t352;
	t338 = t308 * t352;
	t336 = t309 * t351;
	t335 = t363 * t310 * t338 + t353 * t346;
	t330 = t309 * t356 + t359;
	t318 = cos(qJ(3));
	t328 = -t318 * pkin(3) - pkin(4) * t311 - pkin(2) - t347;
	t327 = (t363 * t370 + (pkin(4) + t369) * t354) * t310;
	t326 = -t363 * t317 - t369;
	t325 = t308 * t351 + t315 * t353;
	t324 = -pkin(4) + t326;
	t323 = t324 * t310;
	t322 = t368 * t311 + t312 * t323 - t363 * t340 + t364 * t357;
	t321 = -t368 * t310 + (t324 * t311 - t347) * t312;
	t320 = -t318 * t360 + t321;
	t319 = -pkin(8) - pkin(7);
	t280 = t330 * qJD(1) - t308 * t344 - t311 * t338 - t336;
	t279 = -t308 * t345 + t325 * t311 - t370;
	t278 = t311 * t342 + (t340 + t344) * t309 - t325;
	t277 = t309 * t345 - t311 * t336 - t338 + (t309 * t317 + t311 * t359) * qJD(1);
	t1 = [-t309 * t349 - t363 * t280 - t361 * t279 - t366 * t308 + (-cos(qJ(1)) * pkin(1) + t308 * t319 + t328 * t309) * qJD(1), 0, (-t346 + t362) * t354 + t320 * t309 + t327, t321 * t309 - t346 * t354 + t327, t330 * qJD(6) + t363 * t277 - t361 * t278, -t277; -t308 * t349 - t363 * t278 - t361 * t277 + t366 * t309 + (-sin(qJ(1)) * pkin(1) - t309 * t319 + t328 * t308) * qJD(1), 0, (t323 - t362) * t353 + t320 * t308 + t335, t321 * t308 + t323 * t353 + t335, -(-t308 * t356 + t309 * t315) * qJD(6) + t361 * t280 - t363 * t279, t279; 0, 0, t322 - t348, t322, (-t363 * t315 + t361 * t317) * t357 + (t326 * qJD(5) + t349) * t310, t310 * t351 + t311 * t355;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end