% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(6));
	t151 = t136 * (pkin(8) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(6));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:32
	% EndTime: 2019-10-10 10:24:32
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (84->29), mult. (244->45), div. (0->0), fcn. (222->6), ass. (0->24)
	t175 = sin(pkin(6));
	t194 = t175 * (pkin(8) + r_i_i_C(1));
	t193 = pkin(2) - r_i_i_C(2);
	t191 = r_i_i_C(3) + qJ(3);
	t177 = sin(qJ(2));
	t178 = sin(qJ(1));
	t190 = t178 * t177;
	t179 = cos(qJ(2));
	t189 = t178 * t179;
	t180 = cos(qJ(1));
	t188 = t180 * t177;
	t187 = t180 * t179;
	t186 = qJD(2) * t177;
	t176 = cos(pkin(6));
	t185 = t176 * t190;
	t184 = t176 * t187;
	t183 = qJD(2) * t176 + qJD(1);
	t182 = t176 * t189 + t188;
	t181 = t176 * t188 + t189;
	t170 = -qJD(1) * t185 - t178 * t186 + t183 * t187;
	t169 = t182 * qJD(1) + t181 * qJD(2);
	t168 = t181 * qJD(1) + t182 * qJD(2);
	t167 = -qJD(1) * t184 - qJD(2) * t187 + t183 * t190;
	t1 = [-(-t184 + t190) * qJD(3) - t193 * t170 - t191 * t169 + (-t180 * pkin(1) - t178 * t194) * qJD(1), -(t185 - t187) * qJD(3) - t191 * t168 + t193 * t167, -t167, 0, 0, 0; t182 * qJD(3) - t193 * t168 - t191 * t167 + (-t178 * pkin(1) + t180 * t194) * qJD(1), t181 * qJD(3) - t193 * t169 + t191 * t170, t169, 0, 0, 0; 0, (t177 * qJD(3) + (-t193 * t177 + t191 * t179) * qJD(2)) * t175, t175 * t186, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:32
	% EndTime: 2019-10-10 10:24:32
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (179->58), mult. (534->98), div. (0->0), fcn. (502->8), ass. (0->43)
	t275 = pkin(3) + pkin(8);
	t243 = cos(pkin(6));
	t245 = sin(qJ(2));
	t249 = cos(qJ(1));
	t267 = t249 * t245;
	t246 = sin(qJ(1));
	t248 = cos(qJ(2));
	t268 = t246 * t248;
	t235 = t243 * t268 + t267;
	t251 = t243 * t267 + t268;
	t231 = t235 * qJD(1) + t251 * qJD(2);
	t244 = sin(qJ(4));
	t274 = t231 * t244;
	t247 = cos(qJ(4));
	t273 = t231 * t247;
	t242 = sin(pkin(6));
	t272 = t242 * t246;
	t271 = t242 * t248;
	t270 = t242 * t249;
	t269 = t246 * t245;
	t266 = t249 * t248;
	t265 = qJD(1) * t246;
	t264 = qJD(1) * t249;
	t263 = qJD(2) * t245;
	t259 = t243 * t266;
	t233 = -t259 + t269;
	t262 = qJD(4) * t233;
	t261 = -r_i_i_C(3) - pkin(9) - pkin(2);
	t260 = t243 * t269;
	t258 = t242 * t265;
	t257 = t242 * t264;
	t256 = t242 * t263;
	t255 = qJD(2) * t243 + qJD(1);
	t254 = r_i_i_C(1) * t247 - r_i_i_C(2) * t244;
	t253 = -t244 * r_i_i_C(1) - t247 * r_i_i_C(2);
	t252 = qJ(3) - t253;
	t250 = t254 * qJD(4) + qJD(3);
	t232 = -qJD(1) * t260 - t246 * t263 + t255 * t266;
	t230 = t251 * qJD(1) + t235 * qJD(2);
	t229 = -qJD(1) * t259 - qJD(2) * t266 + t255 * t269;
	t228 = t247 * t257 - t229 * t244 + (t235 * t247 - t244 * t272) * qJD(4);
	t227 = -t244 * t257 - t229 * t247 + (-t235 * t244 - t247 * t272) * qJD(4);
	t1 = [(-t247 * t262 - t274) * r_i_i_C(1) + (t244 * t262 - t273) * r_i_i_C(2) - t231 * qJ(3) - t233 * qJD(3) - pkin(1) * t264 + t261 * t232 + (t253 * t249 * qJD(4) + (-t254 - t275) * t265) * t242, t250 * (-t260 + t266) - t252 * t230 - t261 * t229, -t229, t227 * r_i_i_C(1) - t228 * r_i_i_C(2), 0, 0; t228 * r_i_i_C(1) + t227 * r_i_i_C(2) - t229 * qJ(3) + t235 * qJD(3) + t261 * t230 + (-pkin(1) * t246 + t275 * t270) * qJD(1), t261 * t231 + t252 * t232 + t250 * t251, t231, (-t244 * t258 + t273) * r_i_i_C(1) + (-t247 * t258 - t274) * r_i_i_C(2) + ((-t233 * t244 + t247 * t270) * r_i_i_C(1) + (-t233 * t247 - t244 * t270) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t250 * t245 + (t261 * t245 + t252 * t248) * qJD(2)) * t242, t256, t254 * t256 + ((-t243 * t247 + t244 * t271) * r_i_i_C(1) + (t243 * t244 + t247 * t271) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:32
	% EndTime: 2019-10-10 10:24:32
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (298->82), mult. (727->130), div. (0->0), fcn. (682->10), ass. (0->50)
	t262 = cos(qJ(4));
	t292 = t262 * pkin(4);
	t291 = pkin(8) + pkin(3) + t292;
	t257 = cos(pkin(6));
	t260 = sin(qJ(2));
	t264 = cos(qJ(1));
	t283 = t264 * t260;
	t261 = sin(qJ(1));
	t263 = cos(qJ(2));
	t284 = t261 * t263;
	t243 = t257 * t283 + t284;
	t244 = t257 * t284 + t283;
	t240 = t244 * qJD(1) + t243 * qJD(2);
	t255 = qJ(4) + pkin(11);
	t253 = sin(t255);
	t290 = t240 * t253;
	t254 = cos(t255);
	t289 = t240 * t254;
	t256 = sin(pkin(6));
	t288 = t256 * t261;
	t287 = t256 * t263;
	t286 = t256 * t264;
	t285 = t261 * t260;
	t282 = t264 * t263;
	t281 = qJD(1) * t261;
	t280 = qJD(1) * t264;
	t279 = qJD(2) * t260;
	t278 = qJD(2) * t263;
	t274 = t257 * t282;
	t242 = -t274 + t285;
	t277 = qJD(4) * t242;
	t276 = -r_i_i_C(3) - qJ(5) - pkin(9) - pkin(2);
	t275 = t257 * t285;
	t273 = t256 * t281;
	t272 = t256 * t280;
	t271 = t256 * t279;
	t270 = qJD(2) * t257 + qJD(1);
	t269 = r_i_i_C(1) * t254 - r_i_i_C(2) * t253;
	t259 = sin(qJ(4));
	t268 = -t259 * pkin(4) - t253 * r_i_i_C(1) - t254 * r_i_i_C(2);
	t267 = t269 + t292;
	t266 = qJ(3) - t268;
	t265 = t267 * qJD(4) + qJD(3);
	t245 = -t275 + t282;
	t241 = -qJD(1) * t275 - t261 * t279 + t270 * t282;
	t239 = t243 * qJD(1) + t244 * qJD(2);
	t238 = -qJD(1) * t274 - t264 * t278 + t270 * t285;
	t237 = t254 * t272 - t238 * t253 + (t244 * t254 - t253 * t288) * qJD(4);
	t236 = -t253 * t272 - t238 * t254 + (-t244 * t253 - t254 * t288) * qJD(4);
	t1 = [(-t254 * t277 - t290) * r_i_i_C(1) + (t253 * t277 - t289) * r_i_i_C(2) - t243 * qJD(5) - t240 * qJ(3) - t242 * qJD(3) - pkin(1) * t280 + (-t240 * t259 - t262 * t277) * pkin(4) + t276 * t241 + (t268 * t264 * qJD(4) + (-t269 - t291) * t281) * t256, -t244 * qJD(5) - t276 * t238 - t266 * t239 + t265 * t245, -t238, t236 * r_i_i_C(1) - t237 * r_i_i_C(2) + (-t259 * t272 - t238 * t262 + (-t244 * t259 - t262 * t288) * qJD(4)) * pkin(4), -t239, 0; t237 * r_i_i_C(1) + t236 * r_i_i_C(2) - t238 * qJ(3) + t244 * qJD(3) + t245 * qJD(5) + t276 * t239 + (-pkin(1) * t261 + t291 * t286) * qJD(1) + (-t238 * t259 + (t244 * t262 - t259 * t288) * qJD(4)) * pkin(4), -t242 * qJD(5) + t276 * t240 + t266 * t241 + t265 * t243, t240, (-t253 * t273 + t289) * r_i_i_C(1) + (-t254 * t273 - t290) * r_i_i_C(2) + ((-t242 * t253 + t254 * t286) * r_i_i_C(1) + (-t242 * t254 - t253 * t286) * r_i_i_C(2)) * qJD(4) + (-t259 * t273 + t240 * t262 + (-t242 * t259 + t262 * t286) * qJD(4)) * pkin(4), t241, 0; 0, (qJD(5) * t263 + t265 * t260 + (t276 * t260 + t266 * t263) * qJD(2)) * t256, t271, t267 * t271 + ((t253 * t287 - t254 * t257) * r_i_i_C(1) + (t253 * t257 + t254 * t287) * r_i_i_C(2) + (-t257 * t262 + t259 * t287) * pkin(4)) * qJD(4), t256 * t278, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:34
	% EndTime: 2019-10-10 10:24:35
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (732->117), mult. (1617->183), div. (0->0), fcn. (1604->12), ass. (0->68)
	t398 = cos(pkin(6));
	t406 = cos(qJ(2));
	t407 = cos(qJ(1));
	t435 = t407 * t406;
	t429 = t398 * t435;
	t402 = sin(qJ(2));
	t403 = sin(qJ(1));
	t438 = t403 * t402;
	t381 = -t429 + t438;
	t396 = qJ(4) + pkin(11);
	t394 = sin(t396);
	t395 = cos(t396);
	t397 = sin(pkin(6));
	t439 = t397 * t407;
	t430 = t395 * t439;
	t377 = -t381 * t394 + t430;
	t436 = t407 * t402;
	t437 = t403 * t406;
	t382 = t398 * t436 + t437;
	t400 = sin(qJ(6));
	t404 = cos(qJ(6));
	t451 = -t377 * t400 - t382 * t404;
	t450 = t377 * t404 - t382 * t400;
	t383 = t398 * t437 + t436;
	t370 = t383 * qJD(1) + t382 * qJD(2);
	t434 = qJD(1) * t397;
	t428 = t403 * t434;
	t449 = (qJD(4) * t381 + t428) * t394 - qJD(4) * t430 - t370 * t395;
	t418 = t381 * t395 + t394 * t439;
	t363 = t418 * qJD(4) + t370 * t394 + t395 * t428;
	t421 = t404 * r_i_i_C(1) - t400 * r_i_i_C(2);
	t410 = -t421 * qJD(6) - qJD(5);
	t401 = sin(qJ(4));
	t419 = pkin(5) + t421;
	t448 = pkin(10) + r_i_i_C(3);
	t409 = t401 * pkin(4) + t419 * t394 - t448 * t395 + qJ(3);
	t405 = cos(qJ(4));
	t447 = t405 * pkin(4);
	t446 = -pkin(2) - qJ(5) - pkin(9);
	t445 = pkin(8) + pkin(3) + t447;
	t442 = t397 * t402;
	t441 = t397 * t403;
	t440 = t397 * t406;
	t433 = qJD(2) * t402;
	t432 = qJD(2) * t406;
	t431 = t398 * t438;
	t427 = t407 * t434;
	t426 = t397 * t433;
	t425 = t397 * t432;
	t422 = qJD(2) * t398 + qJD(1);
	t420 = -t400 * r_i_i_C(1) - t404 * r_i_i_C(2);
	t375 = t383 * t394 + t395 * t441;
	t417 = t383 * t395 - t394 * t441;
	t416 = t398 * t394 + t395 * t440;
	t415 = t394 * t440 - t398 * t395;
	t414 = qJD(6) * t420;
	t413 = -t420 - t446;
	t408 = qJD(3) + t394 * t414 + (t448 * t394 + t419 * t395 + t447) * qJD(4);
	t384 = -t431 + t435;
	t372 = t416 * qJD(4) - t394 * t426;
	t371 = -qJD(1) * t431 - t403 * t433 + t422 * t435;
	t369 = t382 * qJD(1) + t383 * qJD(2);
	t368 = -qJD(1) * t429 - t407 * t432 + t422 * t438;
	t366 = t417 * qJD(4) - t368 * t394 + t395 * t427;
	t365 = t375 * qJD(4) + t368 * t395 + t394 * t427;
	t360 = t366 * t404 - t369 * t400 + (-t375 * t400 + t384 * t404) * qJD(6);
	t359 = -t366 * t400 - t369 * t404 + (-t375 * t404 - t384 * t400) * qJD(6);
	t1 = [-t370 * qJ(3) - t381 * qJD(3) - t382 * qJD(5) - t419 * t363 - t413 * t371 - t448 * t449 + (t451 * r_i_i_C(1) - t450 * r_i_i_C(2)) * qJD(6) + (-t407 * pkin(1) - t445 * t441) * qJD(1) + (-t370 * t401 + (-t381 * t405 - t401 * t439) * qJD(4)) * pkin(4), t413 * t368 - t409 * t369 + t410 * t383 + t408 * t384, -t368, t448 * t366 + t417 * t414 - t419 * t365 + (-t401 * t427 - t368 * t405 + (-t383 * t401 - t405 * t441) * qJD(4)) * pkin(4), -t369, t359 * r_i_i_C(1) - t360 * r_i_i_C(2); t366 * pkin(5) + t360 * r_i_i_C(1) + t359 * r_i_i_C(2) - t368 * qJ(3) + t383 * qJD(3) + t384 * qJD(5) + t446 * t369 + t448 * t365 + (-pkin(1) * t403 + t445 * t439) * qJD(1) + (-t368 * t401 + (t383 * t405 - t401 * t441) * qJD(4)) * pkin(4), -t413 * t370 + t409 * t371 + t410 * t381 + t408 * t382, t370, t448 * t363 + t418 * t414 - t419 * t449 + (-t401 * t428 + t370 * t405 + (-t381 * t401 + t405 * t439) * qJD(4)) * pkin(4), t371, (-t363 * t400 + t371 * t404) * r_i_i_C(1) + (-t363 * t404 - t371 * t400) * r_i_i_C(2) + (t450 * r_i_i_C(1) + t451 * r_i_i_C(2)) * qJD(6); 0, ((t409 * qJD(2) - t410) * t406 + (-t413 * qJD(2) + t408) * t402) * t397, t426, -t448 * t372 - t416 * t414 + t419 * (t415 * qJD(4) + t395 * t426) + (t405 * t426 + (-t398 * t405 + t401 * t440) * qJD(4)) * pkin(4), t425, (t372 * t400 + t404 * t425) * r_i_i_C(1) + (t372 * t404 - t400 * t425) * r_i_i_C(2) + ((-t400 * t442 + t404 * t415) * r_i_i_C(1) + (-t400 * t415 - t404 * t442) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end