% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:00
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
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
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
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
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
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
	% StartTime: 2019-10-10 12:00:37
	% EndTime: 2019-10-10 12:00:37
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (179->38), mult. (217->49), div. (0->0), fcn. (148->6), ass. (0->35)
	t189 = qJ(2) + qJ(3);
	t187 = cos(t189);
	t218 = r_i_i_C(3) + qJ(4);
	t202 = t218 * t187;
	t186 = sin(t189);
	t184 = t186 * qJD(4);
	t188 = qJD(2) + qJD(3);
	t190 = sin(qJ(2));
	t217 = pkin(2) * qJD(2);
	t210 = t190 * t217;
	t221 = pkin(3) - r_i_i_C(2);
	t226 = (-t221 * t186 + t202) * t188 + (r_i_i_C(1) + pkin(8) + pkin(7)) * qJD(1) + t184 - t210;
	t191 = sin(qJ(1));
	t214 = t188 * t191;
	t209 = t187 * t214;
	t193 = cos(qJ(1));
	t212 = qJD(1) * t193;
	t225 = t186 * t212 + t209;
	t220 = pkin(2) * t190;
	t216 = t187 * t188;
	t215 = t188 * t186;
	t213 = qJD(1) * t191;
	t211 = qJD(4) * t187;
	t208 = t193 * t216;
	t205 = t186 * t213;
	t207 = pkin(3) * t205 + r_i_i_C(2) * t208 + t193 * t211;
	t203 = t218 * t186;
	t201 = t225 * r_i_i_C(2) + t191 * t211 + t212 * t202;
	t200 = -r_i_i_C(2) * t186 - t202;
	t198 = -t221 * t215 + t218 * t216 + t184;
	t197 = (-pkin(3) * t187 - t203) * t188;
	t192 = cos(qJ(2));
	t196 = qJD(1) * (-t192 * pkin(2) - t221 * t187 - pkin(1) - t203);
	t195 = -t192 * t217 + t197;
	t1 = [-t226 * t191 + t193 * t196, t195 * t193 + (t200 + t220) * t213 + t207, t193 * t197 + t200 * t213 + t207, -t205 + t208, 0, 0; t191 * t196 + t226 * t193, (-pkin(3) * t186 - t220) * t212 + t195 * t191 + t201, -pkin(3) * t209 + (-pkin(3) * t212 - t218 * t214) * t186 + t201, t225, 0, 0; 0, t198 - t210, t198, t215, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:37
	% EndTime: 2019-10-10 12:00:37
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (334->60), mult. (461->87), div. (0->0), fcn. (348->8), ass. (0->55)
	t258 = qJ(2) + qJ(3);
	t256 = cos(t258);
	t259 = sin(qJ(5));
	t262 = cos(qJ(5));
	t310 = r_i_i_C(1) * t259 + r_i_i_C(2) * t262 + qJ(4);
	t314 = t256 * t310;
	t257 = qJD(2) + qJD(3);
	t301 = t256 * t257;
	t250 = qJ(4) * t301;
	t255 = sin(t258);
	t292 = pkin(3) + pkin(9) + r_i_i_C(3);
	t281 = t292 * t257;
	t260 = sin(qJ(2));
	t302 = pkin(2) * qJD(2);
	t291 = t260 * t302;
	t313 = (qJD(4) - t281) * t255 + (pkin(4) + pkin(8) + pkin(7)) * qJD(1) + t250 - t291;
	t293 = qJD(5) * t262;
	t284 = t256 * t293;
	t312 = r_i_i_C(1) * t284 + qJD(4) * t256;
	t279 = qJD(5) * t255 + qJD(1);
	t309 = t262 * t279;
	t308 = t279 * t259;
	t306 = pkin(2) * t260;
	t300 = t257 * t259;
	t299 = t257 * t262;
	t264 = cos(qJ(1));
	t298 = t257 * t264;
	t261 = sin(qJ(1));
	t297 = qJD(1) * t261;
	t296 = qJD(1) * t264;
	t294 = qJD(5) * t259;
	t290 = t256 * t299;
	t289 = t261 * t301;
	t288 = t256 * t298;
	t286 = t255 * t297;
	t285 = t256 * t294;
	t283 = t292 * t255;
	t282 = t292 * t256;
	t278 = -qJD(1) * t255 - qJD(5);
	t277 = t312 * t261 + t296 * t314;
	t276 = t312 * t264 + t292 * t286;
	t275 = t278 * t264;
	t272 = t310 * t255;
	t271 = -r_i_i_C(2) * t294 - t281;
	t263 = cos(qJ(2));
	t270 = qJD(1) * (-t263 * pkin(2) - qJ(4) * t255 - pkin(1) - t282);
	t269 = t278 * t261 + t288;
	t268 = t256 * r_i_i_C(1) * t300 + r_i_i_C(2) * t290 + t250 + (r_i_i_C(1) * t293 + qJD(4) + t271) * t255;
	t267 = -r_i_i_C(2) * t285 + (-t282 - t272) * t257;
	t266 = -t263 * t302 + t267;
	t238 = t269 * t259 + t264 * t309;
	t237 = t269 * t262 - t264 * t308;
	t236 = -t261 * t309 + (t275 - t289) * t259;
	t235 = t262 * t275 + (-t290 + t308) * t261;
	t1 = [t236 * r_i_i_C(1) + t235 * r_i_i_C(2) - t313 * t261 + t264 * t270, (t306 - t314) * t297 + t266 * t264 + t276, -t272 * t298 + (t271 * t264 - t297 * t310) * t256 + t276, -t286 + t288, t237 * r_i_i_C(1) - t238 * r_i_i_C(2), 0; t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t261 * t270 + t313 * t264, (-t283 - t306) * t296 + t266 * t261 + t277, t267 * t261 - t283 * t296 + t277, t255 * t296 + t289, -t235 * r_i_i_C(1) + t236 * r_i_i_C(2), 0; 0, t268 - t291, t268, t257 * t255, (-t255 * t300 + t284) * r_i_i_C(2) + (t255 * t299 + t285) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:37
	% EndTime: 2019-10-10 12:00:38
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (610->80), mult. (663->111), div. (0->0), fcn. (506->10), ass. (0->68)
	t315 = qJ(2) + qJ(3);
	t311 = cos(t315);
	t314 = qJ(5) + qJ(6);
	t308 = sin(t314);
	t310 = cos(t314);
	t316 = sin(qJ(5));
	t342 = pkin(5) * t316 + qJ(4);
	t379 = r_i_i_C(1) * t308 + r_i_i_C(2) * t310 + t342;
	t386 = t311 * t379;
	t319 = cos(qJ(5));
	t368 = t319 * pkin(5);
	t382 = qJD(5) * t368 + qJD(4);
	t385 = pkin(3) + r_i_i_C(3);
	t309 = sin(t315);
	t321 = cos(qJ(1));
	t356 = qJD(1) * t321;
	t343 = t309 * t356;
	t318 = sin(qJ(1));
	t313 = qJD(2) + qJD(3);
	t365 = t311 * t313;
	t347 = t318 * t365;
	t384 = t343 + t347;
	t312 = qJD(5) + qJD(6);
	t366 = t311 * t312;
	t348 = t310 * t366;
	t381 = r_i_i_C(1) * t348 + t382 * t311;
	t317 = sin(qJ(2));
	t367 = pkin(2) * qJD(2);
	t350 = t317 * t367;
	t322 = -pkin(10) - pkin(9);
	t353 = -t322 + t385;
	t380 = (t353 * t313 - t382) * t309 - (pkin(4) + t368 + pkin(8) + pkin(7)) * qJD(1) - t342 * t365 + t350;
	t340 = t309 * t312 + qJD(1);
	t378 = t318 * t340;
	t377 = t321 * t340;
	t373 = pkin(2) * t317;
	t370 = r_i_i_C(2) * t308;
	t364 = t313 * t309;
	t363 = t313 * t319;
	t362 = t313 * t321;
	t358 = qJD(1) * t309;
	t339 = -t312 - t358;
	t329 = t339 * t321 - t347;
	t276 = t308 * t378 + t329 * t310;
	t277 = t329 * t308 - t310 * t378;
	t361 = -t276 * r_i_i_C(1) + t277 * r_i_i_C(2);
	t346 = t311 * t362;
	t328 = t339 * t318 + t346;
	t278 = -t308 * t377 + t328 * t310;
	t279 = t328 * t308 + t310 * t377;
	t360 = t278 * r_i_i_C(1) - t279 * r_i_i_C(2);
	t357 = qJD(1) * t318;
	t354 = qJD(5) * t316;
	t352 = r_i_i_C(1) * t309 * t310;
	t349 = t308 * t366;
	t344 = t309 * t357;
	t337 = qJD(5) + t358;
	t335 = (-qJD(5) * t309 - qJD(1)) * t316;
	t334 = t381 * t321 + t322 * t346 + t385 * t344;
	t333 = -t312 * t370 - t313 * t385;
	t332 = r_i_i_C(1) * t349 + r_i_i_C(2) * t348 + t313 * t352 - t364 * t370;
	t330 = t381 * t318 + t384 * t322 + t356 * t386;
	t320 = cos(qJ(2));
	t327 = -pkin(5) * t354 + (-t320 * pkin(2) - t342 * t309 - t353 * t311 - pkin(1)) * qJD(1);
	t326 = t312 * t352 + t322 * t364 + (t333 + t382) * t309 + t379 * t365;
	t325 = -r_i_i_C(2) * t349 + (-t309 * t379 - t311 * t385) * t313;
	t324 = -t320 * t367 + t325;
	t1 = [t277 * r_i_i_C(1) + t276 * r_i_i_C(2) + t380 * t318 + t327 * t321, (-t309 * t322 + t373 - t386) * t357 + t324 * t321 + t334, (-t322 * t357 - t362 * t379) * t309 + (t333 * t321 - t357 * t379) * t311 + t334, -t344 + t346, (t321 * t335 + (-t337 * t318 + t346) * t319) * pkin(5) + t360, t360; t279 * r_i_i_C(1) + t278 * r_i_i_C(2) + t327 * t318 - t380 * t321, (-t309 * t385 - t373) * t356 + t324 * t318 + t330, t325 * t318 - t343 * t385 + t330, t384, (t337 * t321 * t319 + (t311 * t363 + t335) * t318) * pkin(5) + t361, t361; 0, t326 - t350, t326, t364, (t309 * t363 + t311 * t354) * pkin(5) + t332, t332;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end