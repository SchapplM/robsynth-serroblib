% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRP5
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
% Datum: 2019-10-10 01:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:23
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
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(10)) + r_i_i_C(2) * sin(pkin(10)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:23
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(7) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(10) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(10)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t28, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:23
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (124->30), mult. (120->40), div. (0->0), fcn. (79->7), ass. (0->28)
	t44 = qJD(3) + qJD(4);
	t43 = pkin(10) + qJ(3);
	t41 = qJ(4) + t43;
	t38 = cos(t41);
	t61 = r_i_i_C(2) * t38;
	t37 = sin(t41);
	t63 = r_i_i_C(1) * t37;
	t51 = t61 + t63;
	t49 = t51 * t44;
	t39 = sin(t43);
	t64 = pkin(3) * t39;
	t65 = qJD(3) * t64 + t49;
	t62 = r_i_i_C(2) * t37;
	t60 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(2);
	t59 = t38 * t44;
	t45 = sin(qJ(1));
	t58 = qJD(1) * t45;
	t46 = cos(qJ(1));
	t57 = qJD(1) * t46;
	t40 = cos(t43);
	t56 = qJD(3) * t40;
	t55 = r_i_i_C(1) * t59;
	t54 = t44 * t62;
	t52 = qJD(1) * t61;
	t50 = -r_i_i_C(1) * t38 - pkin(3) * t40 - cos(pkin(10)) * pkin(2) - pkin(1) + t62;
	t48 = t45 * t52 + t58 * t63 + (t54 - t55) * t46;
	t32 = t45 * t54;
	t1 = [t46 * qJD(2) + t65 * t45 + (-t60 * t45 + t50 * t46) * qJD(1), t57, (t39 * t58 - t46 * t56) * pkin(3) + t48, t48, 0, 0; t45 * qJD(2) - t65 * t46 + (t50 * t45 + t60 * t46) * qJD(1), t58, t32 + (-pkin(3) * t56 - t55) * t45 + (-t51 - t64) * t57, -t46 * t52 + t32 + (-t37 * t57 - t45 * t59) * r_i_i_C(1), 0, 0; 0, 0, -t65, -t49, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:24
	% EndTime: 2019-10-10 01:52:25
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (387->62), mult. (400->98), div. (0->0), fcn. (303->9), ass. (0->55)
	t262 = pkin(10) + qJ(3);
	t260 = qJ(4) + t262;
	t256 = sin(t260);
	t266 = cos(qJ(5));
	t312 = r_i_i_C(1) * t266 + pkin(4);
	t315 = t256 * t312;
	t264 = sin(qJ(5));
	t295 = qJD(5) * t266;
	t257 = cos(t260);
	t263 = qJD(3) + qJD(4);
	t303 = t257 * t263;
	t314 = t256 * t295 + t264 * t303;
	t309 = pkin(9) + r_i_i_C(3);
	t290 = t309 * t257;
	t258 = sin(t262);
	t304 = pkin(3) * qJD(3);
	t293 = t258 * t304;
	t307 = pkin(4) * t256;
	t313 = -t293 + (t290 - t307) * t263;
	t296 = qJD(5) * t264;
	t283 = t256 * t296;
	t310 = r_i_i_C(1) * t283 + t314 * r_i_i_C(2);
	t308 = pkin(3) * t258;
	t305 = r_i_i_C(2) * t264;
	t265 = sin(qJ(1));
	t302 = t263 * t265;
	t301 = t263 * t266;
	t267 = cos(qJ(1));
	t300 = t263 * t267;
	t299 = t266 * t267;
	t298 = qJD(1) * t265;
	t297 = qJD(1) * t267;
	t294 = t256 * t305;
	t292 = qJD(1) * t305;
	t291 = t309 * t256;
	t289 = t309 * t265;
	t288 = t256 * t301;
	t278 = qJD(5) * t257 - qJD(1);
	t277 = qJD(1) * t257 - qJD(5);
	t276 = t312 * t257;
	t275 = t312 * t267;
	t274 = t310 * t267 + t298 * t315;
	t273 = t278 * t264;
	t272 = t267 * t256 * t292 + t310 * t265 + t297 * t290;
	t259 = cos(t262);
	t271 = -pkin(4) * t257 - pkin(3) * t259 - cos(pkin(10)) * pkin(2) - pkin(1) - t291;
	t270 = t256 * t300 + t277 * t265;
	t269 = -t259 * t304 + (-t276 - t291) * t263;
	t268 = -t257 * r_i_i_C(2) * t295 + (-t257 * t296 - t288) * r_i_i_C(1) + t309 * t303 + (-t307 + t294) * t263;
	t261 = -pkin(8) - pkin(7) - qJ(2);
	t240 = -t277 * t299 + (t273 + t288) * t265;
	t239 = t278 * t266 * t265 + (-t256 * t302 + t277 * t267) * t264;
	t238 = t270 * t266 + t267 * t273;
	t237 = t270 * t264 - t278 * t299;
	t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) + t267 * qJD(2) - t313 * t265 + (t261 * t265 + t271 * t267) * qJD(1), t297, (-t290 - t294 + t308) * t298 + t269 * t267 + t274, (-t265 * t292 - t309 * t300) * t256 + (-qJD(1) * t289 - t263 * t275) * t257 + t274, t237 * r_i_i_C(1) + t238 * r_i_i_C(2), 0; -t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t265 * qJD(2) + t313 * t267 + (-t261 * t267 + t271 * t265) * qJD(1), t298, (-t308 - t315) * t297 + t269 * t265 + t272, -t276 * t302 + (-qJD(1) * t275 - t263 * t289) * t256 + t272, -t239 * r_i_i_C(1) + t240 * r_i_i_C(2), 0; 0, 0, t268 - t293, t268, (-t257 * t301 + t283) * r_i_i_C(2) - t314 * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:25
	% EndTime: 2019-10-10 01:52:25
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (670->74), mult. (708->100), div. (0->0), fcn. (571->9), ass. (0->58)
	t321 = sin(qJ(5));
	t323 = cos(qJ(5));
	t324 = cos(qJ(1));
	t360 = qJD(5) * t324;
	t322 = sin(qJ(1));
	t363 = qJD(1) * t322;
	t332 = t321 * t360 + t323 * t363;
	t373 = pkin(5) + r_i_i_C(1);
	t370 = r_i_i_C(3) + qJ(6);
	t377 = t370 * t321;
	t319 = pkin(10) + qJ(3);
	t317 = qJ(4) + t319;
	t314 = cos(t317);
	t372 = pkin(9) + r_i_i_C(2);
	t354 = t372 * t314;
	t320 = qJD(3) + qJD(4);
	t353 = t372 * t320;
	t315 = sin(t319);
	t369 = pkin(3) * qJD(3);
	t357 = t315 * t369;
	t359 = qJD(6) * t321;
	t313 = sin(t317);
	t368 = t313 * t320;
	t376 = -pkin(4) * t368 + (t353 + t359) * t314 - t357;
	t361 = qJD(5) * t323;
	t330 = -t370 * t361 - t359;
	t371 = pkin(3) * t315;
	t367 = t314 * t320;
	t366 = t314 * t322;
	t365 = t322 * t321;
	t364 = t324 * t323;
	t362 = qJD(1) * t324;
	t358 = t323 * qJD(6);
	t356 = t373 * t321;
	t355 = t372 * t313;
	t352 = t322 * t368;
	t351 = t324 * t368;
	t346 = qJD(5) * t365;
	t344 = t323 * t360;
	t343 = qJD(2) - t358;
	t342 = t373 * t313 * t346 + t362 * t354;
	t337 = t314 * t364 + t365;
	t316 = cos(t319);
	t335 = -pkin(4) * t314 - pkin(3) * t316 - cos(pkin(10)) * pkin(2) - pkin(1) - t355;
	t334 = (t373 * t332 + (pkin(4) + t377) * t363) * t313;
	t333 = -t373 * t323 - t377;
	t331 = t321 * t362 + t322 * t361;
	t329 = -pkin(4) + t333;
	t328 = t329 * t320;
	t327 = t313 * t328 + t372 * t367 + (-qJD(5) * t356 - t330) * t314;
	t326 = t330 * t313 + (t329 * t314 - t355) * t320;
	t325 = -t316 * t369 + t326;
	t318 = -pkin(8) - pkin(7) - qJ(2);
	t285 = t337 * qJD(1) - t314 * t346 - t323 * t352 - t344;
	t284 = t331 * t314 - t321 * t352 - t332;
	t283 = t332 * t314 + t323 * t351 - t331;
	t282 = t321 * t351 - t314 * t344 - t346 + (t314 * t365 + t364) * qJD(1);
	t1 = [t343 * t324 - t373 * t285 - t370 * t284 - t376 * t322 + (t322 * t318 + t335 * t324) * qJD(1), t362, (-t354 + t371) * t363 + t325 * t324 + t334, t326 * t324 - t354 * t363 + t334, t337 * qJD(6) + t373 * t282 - t370 * t283, -t282; t343 * t322 - t373 * t283 - t370 * t282 + t376 * t324 + (-t324 * t318 + t335 * t322) * qJD(1), t363, (t329 * t313 - t371) * t362 + t325 * t322 + t342, t328 * t366 + ((-t353 + t330) * t322 + t329 * t362) * t313 + t342, -(t324 * t321 - t323 * t366) * qJD(6) + t370 * t285 - t373 * t284, t284; 0, 0, t327 - t357, t327, (t370 * t323 - t356) * t367 + (t333 * qJD(5) + t358) * t313, t313 * t361 + t321 * t367;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end