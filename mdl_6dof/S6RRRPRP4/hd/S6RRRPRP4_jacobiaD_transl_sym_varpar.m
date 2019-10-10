% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:40
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
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
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
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
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
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
	% StartTime: 2019-10-10 11:40:27
	% EndTime: 2019-10-10 11:40:27
	% DurationCPUTime: 0.23s
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
	t226 = (-t186 * t221 + t202) * t188 + (r_i_i_C(1) + pkin(8) + pkin(7)) * qJD(1) + t184 - t210;
	t191 = sin(qJ(1));
	t214 = t188 * t191;
	t209 = t187 * t214;
	t193 = cos(qJ(1));
	t212 = qJD(1) * t193;
	t224 = t186 * t212 + t209;
	t220 = pkin(2) * t190;
	t216 = t187 * t188;
	t215 = t188 * t186;
	t213 = qJD(1) * t191;
	t211 = qJD(4) * t187;
	t208 = t193 * t216;
	t205 = t186 * t213;
	t207 = pkin(3) * t205 + r_i_i_C(2) * t208 + t193 * t211;
	t203 = t218 * t186;
	t201 = t224 * r_i_i_C(2) + t191 * t211 + t212 * t202;
	t200 = -r_i_i_C(2) * t186 - t202;
	t198 = -t221 * t215 + t218 * t216 + t184;
	t197 = (-pkin(3) * t187 - t203) * t188;
	t192 = cos(qJ(2));
	t196 = qJD(1) * (-t192 * pkin(2) - t187 * t221 - pkin(1) - t203);
	t195 = -t192 * t217 + t197;
	t1 = [-t226 * t191 + t193 * t196, t195 * t193 + (t200 + t220) * t213 + t207, t193 * t197 + t200 * t213 + t207, -t205 + t208, 0, 0; t191 * t196 + t226 * t193, (-pkin(3) * t186 - t220) * t212 + t195 * t191 + t201, -pkin(3) * t209 + (-pkin(3) * t212 - t214 * t218) * t186 + t201, t224, 0, 0; 0, t198 - t210, t198, t215, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:27
	% EndTime: 2019-10-10 11:40:27
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
	% StartTime: 2019-10-10 11:40:28
	% EndTime: 2019-10-10 11:40:28
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (538->80), mult. (769->116), div. (0->0), fcn. (616->8), ass. (0->65)
	t382 = pkin(5) + r_i_i_C(1);
	t379 = r_i_i_C(3) + qJ(6);
	t362 = pkin(3) + pkin(9) + r_i_i_C(2);
	t315 = qJ(2) + qJ(3);
	t313 = cos(t315);
	t314 = qJD(2) + qJD(3);
	t377 = t313 * t314;
	t306 = qJ(4) * t377;
	t312 = sin(t315);
	t319 = cos(qJ(5));
	t364 = qJD(6) * t319;
	t331 = -t362 * t314 - t364;
	t317 = sin(qJ(2));
	t378 = pkin(2) * qJD(2);
	t361 = t317 * t378;
	t386 = (qJD(4) + t331) * t312 + (pkin(4) + pkin(8) + pkin(7)) * qJD(1) + t306 - t361;
	t316 = sin(qJ(5));
	t385 = t382 * t316;
	t381 = pkin(2) * t317;
	t318 = sin(qJ(1));
	t376 = t313 * t318;
	t375 = t314 * t312;
	t321 = cos(qJ(1));
	t374 = t314 * t321;
	t373 = t316 * t318;
	t372 = t316 * t321;
	t371 = t318 * t319;
	t370 = t321 * t319;
	t369 = qJD(1) * t318;
	t368 = qJD(1) * t321;
	t367 = qJD(4) * t313;
	t366 = qJD(5) * t316;
	t365 = qJD(5) * t319;
	t363 = t316 * qJD(6);
	t359 = t314 * t376;
	t358 = t313 * t374;
	t357 = t314 * t371;
	t356 = t314 * t370;
	t355 = t312 * t369;
	t354 = t312 * t368;
	t353 = t313 * t368;
	t350 = t313 * t366;
	t349 = qJD(5) * t376;
	t348 = t321 * t365;
	t347 = t379 * t319;
	t346 = t362 * t313;
	t345 = qJD(5) * t312 + qJD(1);
	t344 = qJD(1) * t312 + qJD(5);
	t335 = -qJ(4) - t385;
	t334 = t345 * t319;
	t333 = t335 * t312;
	t332 = t312 * t372 + t371;
	t329 = -t362 * t312 - t313 * t347;
	t328 = qJ(4) * t353 + t318 * t367 + t382 * (t316 * t353 + t319 * t349) + t379 * (t312 * t357 + t316 * t349);
	t327 = t382 * t313 * t348 + t321 * t367 + t362 * t355 + t379 * (t313 * t319 * t369 + t312 * t356 + t321 * t350);
	t320 = cos(qJ(2));
	t326 = t363 + (-pkin(2) * t320 - qJ(4) * t312 - pkin(1) - t346) * qJD(1);
	t325 = -t313 * t364 + (-t346 + t333) * t314;
	t324 = -t320 * t378 + t325;
	t323 = t329 * t314 + t306 + t382 * (t312 * t365 + t316 * t377) + (t379 * t366 + qJD(4) - t364) * t312;
	t276 = t321 * t334 + (-t344 * t318 + t358) * t316;
	t275 = -t313 * t356 + t332 * qJD(5) + (t312 * t371 + t372) * qJD(1);
	t274 = t318 * t334 + (t344 * t321 + t359) * t316;
	t273 = -t313 * t357 - t319 * t354 + t345 * t373 - t348;
	t1 = [-t379 * t273 - t382 * t274 - t386 * t318 + t326 * t321, t327 + t324 * t321 + (t335 * t313 + t381) * t369, t327 + (t331 * t321 + t335 * t369) * t313 + t333 * t374, -t355 + t358, t332 * qJD(6) - t382 * t275 + t379 * t276, t275; t379 * t275 + t382 * t276 + t326 * t318 + t386 * t321, (t329 - t381) * t368 + t324 * t318 + t328, t325 * t318 + t329 * t368 + t328, t354 + t359, -(-t312 * t373 + t370) * qJD(6) + t379 * t274 - t382 * t273, t273; 0, t323 - t361, t323, t375, (t379 * t316 + t382 * t319) * t375 + (-t363 + (-t347 + t385) * qJD(5)) * t313, -t319 * t375 - t350;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end