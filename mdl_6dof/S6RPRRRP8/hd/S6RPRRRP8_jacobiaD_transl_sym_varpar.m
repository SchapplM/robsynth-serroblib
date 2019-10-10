% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:55
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
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
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (23->16), mult. (72->29), div. (0->0), fcn. (46->4), ass. (0->13)
	t17 = sin(qJ(3));
	t19 = cos(qJ(3));
	t29 = (r_i_i_C(1) * t19 - r_i_i_C(2) * t17) * qJD(3);
	t18 = sin(qJ(1));
	t28 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t27 = qJD(1) * t20;
	t26 = qJD(3) * t18;
	t25 = qJD(3) * t20;
	t24 = -pkin(1) - pkin(7) - r_i_i_C(3);
	t22 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19 + qJ(2);
	t21 = qJD(2) + t29;
	t1 = [t21 * t20 + (-t18 * t22 + t20 * t24) * qJD(1), t27, (-t17 * t27 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 + t19 * t27) * r_i_i_C(1), 0, 0, 0; t21 * t18 + (t18 * t24 + t20 * t22) * qJD(1), t28, (-t17 * t28 + t19 * t25) * r_i_i_C(2) + (t17 * t25 + t19 * t28) * r_i_i_C(1), 0, 0, 0; 0, 0, -t29, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (85->27), mult. (126->41), div. (0->0), fcn. (81->6), ass. (0->27)
	t43 = cos(qJ(3));
	t62 = pkin(3) * t43;
	t40 = qJ(3) + qJ(4);
	t38 = cos(t40);
	t61 = r_i_i_C(1) * t38;
	t37 = sin(t40);
	t60 = r_i_i_C(2) * t37;
	t39 = qJD(3) + qJD(4);
	t59 = t37 * t39;
	t58 = t38 * t39;
	t42 = sin(qJ(1));
	t57 = qJD(1) * t42;
	t44 = cos(qJ(1));
	t56 = qJD(1) * t44;
	t41 = sin(qJ(3));
	t55 = qJD(3) * t41;
	t54 = -pkin(1) - r_i_i_C(3) - pkin(8) - pkin(7);
	t53 = r_i_i_C(1) * t59;
	t52 = qJD(1) * t61;
	t51 = qJD(3) * t62;
	t50 = -r_i_i_C(1) * t58 + r_i_i_C(2) * t59;
	t49 = -r_i_i_C(1) * t37 - r_i_i_C(2) * t38;
	t48 = t42 * t52 - t57 * t60 + (t58 * r_i_i_C(2) + t53) * t44;
	t47 = pkin(3) * t41 + qJ(2) - t49;
	t46 = t51 + qJD(2) + (-t60 + t61) * t39;
	t35 = t44 * t52;
	t1 = [t46 * t44 + (-t47 * t42 + t54 * t44) * qJD(1), t56, t35 + (-t60 + t62) * t56 + (-pkin(3) * t55 + t49 * t39) * t42, -t42 * t53 + t35 + (-t37 * t56 - t42 * t58) * r_i_i_C(2), 0, 0; t46 * t42 + (t54 * t42 + t47 * t44) * qJD(1), t57, (t43 * t57 + t44 * t55) * pkin(3) + t48, t48, 0, 0; 0, 0, t50 - t51, t50, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:52
	% EndTime: 2019-10-10 01:55:53
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (274->62), mult. (406->96), div. (0->0), fcn. (305->8), ass. (0->56)
	t302 = pkin(9) + r_i_i_C(3);
	t258 = cos(qJ(5));
	t299 = r_i_i_C(1) * t258;
	t308 = pkin(4) + t299;
	t256 = sin(qJ(3));
	t254 = qJ(3) + qJ(4);
	t252 = cos(t254);
	t285 = t302 * t252;
	t251 = sin(t254);
	t301 = pkin(4) * t251;
	t306 = -pkin(3) * t256 - qJ(2) + t285 - t301;
	t260 = cos(qJ(1));
	t271 = qJD(1) * t251 + qJD(5);
	t305 = t271 * t260;
	t304 = t302 * t251;
	t257 = sin(qJ(1));
	t253 = qJD(3) + qJD(4);
	t293 = t253 * t260;
	t303 = -t252 * t293 + t271 * t257;
	t300 = pkin(4) * t252;
	t255 = sin(qJ(5));
	t298 = r_i_i_C(2) * t255;
	t297 = -pkin(1) - pkin(8) - pkin(7);
	t296 = pkin(3) * qJD(3);
	t295 = t251 * t255;
	t294 = t253 * t258;
	t292 = qJD(1) * t257;
	t291 = qJD(1) * t260;
	t290 = qJD(5) * t255;
	t289 = qJD(5) * t258;
	t288 = t252 * t298;
	t287 = t256 * t296;
	t259 = cos(qJ(3));
	t286 = t259 * t296;
	t284 = t253 * t295;
	t282 = t252 * t253 * t257;
	t278 = t252 * t291;
	t277 = t252 * t290;
	t276 = t252 * t289;
	t275 = r_i_i_C(2) * t284;
	t274 = qJD(1) * t252 * t299;
	t273 = r_i_i_C(2) * t276;
	t272 = -qJD(5) * t251 - qJD(1);
	t270 = t272 * t260;
	t268 = pkin(4) * t278 + t257 * t275 + t260 * t274 + t302 * (t251 * t291 + t282);
	t267 = qJD(1) * (pkin(3) * t259 - t288);
	t266 = t251 * t294 + t277;
	t265 = t257 * t274 + t292 * t300 + (r_i_i_C(1) * t277 + t273) * t260 + (t302 * t292 + t308 * t293) * t251;
	t264 = t286 + qJD(2) + (t300 + t304) * t253;
	t263 = (-t308 * t252 + t288 - t304) * t253 + (r_i_i_C(1) * t290 + r_i_i_C(2) * t289) * t251;
	t262 = -t266 * r_i_i_C(1) - t253 * t301 - t273;
	t232 = t258 * t305 + (t252 * t294 + t272 * t255) * t257;
	t231 = t272 * t258 * t257 + (-t282 - t305) * t255;
	t230 = t255 * t270 - t303 * t258;
	t229 = t303 * t255 + t258 * t270;
	t1 = [t230 * r_i_i_C(1) + t229 * r_i_i_C(2) + t264 * t260 + (t306 * t257 + t297 * t260) * qJD(1), t291, t260 * t267 + (t262 - t287) * t257 + t268, t262 * t257 - t278 * t298 + t268, t231 * r_i_i_C(1) - t232 * r_i_i_C(2), 0; t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t264 * t257 + (t297 * t257 - t306 * t260) * qJD(1), t292, t257 * t267 + (t287 + (-r_i_i_C(2) * t295 - t285) * t253) * t260 + t265, -t260 * t275 + (-t292 * t298 - t302 * t293) * t252 + t265, -t229 * r_i_i_C(1) + t230 * r_i_i_C(2), 0; 0, 0, t263 - t286, t263, t266 * r_i_i_C(2) + (-t276 + t284) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:53
	% EndTime: 2019-10-10 01:55:53
	% DurationCPUTime: 0.62s
	% Computational Cost: add. (478->69), mult. (714->101), div. (0->0), fcn. (573->8), ass. (0->60)
	t306 = sin(qJ(5));
	t366 = -r_i_i_C(1) - pkin(5);
	t381 = t366 * t306;
	t380 = qJD(5) * t381;
	t309 = cos(qJ(5));
	t364 = r_i_i_C(3) + qJ(6);
	t321 = -t364 * t306 + t366 * t309;
	t319 = -pkin(4) + t321;
	t305 = qJ(3) + qJ(4);
	t303 = cos(t305);
	t379 = t303 * t380;
	t304 = qJD(3) + qJD(4);
	t352 = qJD(6) * t306;
	t367 = pkin(9) + r_i_i_C(2);
	t378 = t367 * t304 + t352;
	t335 = t364 * t309;
	t376 = -qJD(5) * t335 - t378;
	t374 = t303 * t304;
	t302 = sin(t305);
	t307 = sin(qJ(3));
	t368 = -pkin(3) * t307 - pkin(4) * t302 + t367 * t303 - qJ(2);
	t365 = -pkin(1) - pkin(8) - pkin(7);
	t363 = pkin(3) * qJD(3);
	t362 = t302 * t304;
	t308 = sin(qJ(1));
	t361 = t303 * t308;
	t360 = t303 * t309;
	t311 = cos(qJ(1));
	t359 = t304 * t311;
	t358 = t308 * t306;
	t357 = t308 * t309;
	t356 = t311 * t306;
	t355 = qJD(1) * t308;
	t354 = qJD(1) * t311;
	t351 = t309 * qJD(6);
	t310 = cos(qJ(3));
	t350 = pkin(3) * qJD(1) * t310;
	t349 = t307 * t363;
	t348 = t310 * t363;
	t347 = t304 * t361;
	t346 = t311 * t302 * t309;
	t344 = t303 * t359;
	t340 = t303 * t354;
	t339 = t309 * t354;
	t336 = qJD(5) * t360;
	t334 = qJD(5) * t302 + qJD(1);
	t333 = qJD(1) * t302 + qJD(5);
	t324 = t333 * t311;
	t322 = t302 * t357 + t356;
	t318 = -t366 * t303 * t339 + pkin(4) * t340 + t352 * t361 + t364 * (t306 * t340 + t308 * t336) + t367 * (t302 * t354 + t347);
	t317 = t376 * t303;
	t316 = t367 * t302 * t355 - t311 * t379 - t319 * (t302 * t359 + t303 * t355);
	t315 = pkin(4) * t374 + t378 * t302 + qJD(2) + t348;
	t314 = t319 * t362 + t379;
	t313 = t319 * t374 + (t376 - t380) * t302;
	t269 = t309 * t324 + (t304 * t360 - t334 * t306) * t308;
	t268 = t334 * t357 + (t324 + t347) * t306;
	t267 = -t309 * t344 + (t302 * t356 + t357) * qJD(5) + t322 * qJD(1);
	t266 = -qJD(5) * t346 - t306 * t344 + t333 * t358 - t339;
	t1 = [t308 * t351 + t366 * t267 - t364 * t266 + t315 * t311 + (t368 * t308 + t365 * t311) * qJD(1), t354, (t314 - t349) * t308 + t311 * t350 + t318, t314 * t308 + t318, t322 * qJD(6) + t366 * t268 + t364 * t269, t268; -t311 * t351 - t366 * t269 + t364 * t268 + t315 * t308 + (t365 * t308 - t368 * t311) * qJD(1), t355, t308 * t350 + (t317 + t349) * t311 + t316, t311 * t317 + t316, -(t346 - t358) * qJD(6) + t364 * t267 + t366 * t266, t266; 0, 0, t313 - t348, t313, (-t335 - t381) * t362 + (t321 * qJD(5) + t351) * t303, -t306 * t362 + t336;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end