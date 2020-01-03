% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR10
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:28
	% EndTime: 2019-12-31 21:31:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:28
	% EndTime: 2019-12-31 21:31:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:29
	% EndTime: 2019-12-31 21:31:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(5));
	t151 = t136 * (pkin(7) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(5));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:29
	% EndTime: 2019-12-31 21:31:29
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (145->51), mult. (441->91), div. (0->0), fcn. (412->8), ass. (0->38)
	t242 = cos(pkin(5));
	t244 = sin(qJ(2));
	t245 = sin(qJ(1));
	t264 = t245 * t244;
	t258 = t242 * t264;
	t247 = cos(qJ(2));
	t248 = cos(qJ(1));
	t261 = t248 * t247;
	t233 = -qJD(1) * t258 - qJD(2) * t264 + (qJD(2) * t242 + qJD(1)) * t261;
	t241 = sin(pkin(5));
	t265 = t241 * t248;
	t256 = qJD(3) * t265;
	t269 = t233 - t256;
	t268 = -r_i_i_C(3) - pkin(8);
	t267 = t241 * t245;
	t246 = cos(qJ(3));
	t266 = t241 * t246;
	t263 = t245 * t247;
	t262 = t248 * t244;
	t260 = qJD(1) * t241;
	t235 = t242 * t262 + t263;
	t259 = qJD(3) * t235;
	t257 = t248 * t260;
	t243 = sin(qJ(3));
	t255 = -r_i_i_C(1) * t243 - r_i_i_C(2) * t246;
	t254 = t246 * r_i_i_C(1) - t243 * r_i_i_C(2) + pkin(2);
	t253 = t242 * t261 - t264;
	t252 = t242 * t263 + t262;
	t251 = t258 - t261;
	t250 = qJD(3) * t255;
	t249 = t245 * t260 - t259;
	t238 = t246 * t256;
	t232 = t252 * qJD(1) + t235 * qJD(2);
	t231 = t235 * qJD(1) + t252 * qJD(2);
	t230 = -t253 * qJD(1) + t251 * qJD(2);
	t229 = t243 * t257 - t231 * t246 + (t243 * t251 + t245 * t266) * qJD(3);
	t228 = t246 * t257 + t231 * t243 + (-t243 * t267 + t246 * t251) * qJD(3);
	t1 = [(-t233 * t246 + t243 * t259 + t238) * r_i_i_C(1) + (t269 * t243 + t246 * t259) * r_i_i_C(2) - t233 * pkin(2) + t268 * t232 + (-t248 * pkin(1) + (-pkin(7) + t255) * t267) * qJD(1), t254 * t230 + t268 * t231 - t252 * t250, t228 * r_i_i_C(1) - t229 * r_i_i_C(2), 0, 0; -t231 * pkin(2) + t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t268 * t230 + (-pkin(1) * t245 + pkin(7) * t265) * qJD(1), -t254 * t232 - t268 * t233 + t253 * t250, t238 * r_i_i_C(2) + (t249 * r_i_i_C(1) - t233 * r_i_i_C(2)) * t246 + (-t269 * r_i_i_C(1) - t249 * r_i_i_C(2)) * t243, 0, 0; 0, (t247 * t250 + (-t254 * t244 - t268 * t247) * qJD(2)) * t241, t255 * t247 * t241 * qJD(2) + ((-t242 * t243 - t244 * t266) * r_i_i_C(1) + (t241 * t243 * t244 - t242 * t246) * r_i_i_C(2)) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:29
	% EndTime: 2019-12-31 21:31:30
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (252->77), mult. (593->125), div. (0->0), fcn. (554->10), ass. (0->52)
	t253 = qJ(3) + pkin(10);
	t251 = sin(t253);
	t252 = cos(t253);
	t257 = sin(qJ(3));
	t291 = pkin(3) * t257;
	t264 = r_i_i_C(1) * t251 + r_i_i_C(2) * t252 + t291;
	t263 = qJD(3) * t264;
	t290 = t251 * r_i_i_C(2);
	t260 = cos(qJ(3));
	t289 = t260 * pkin(3);
	t288 = r_i_i_C(3) + qJ(4) + pkin(8);
	t254 = sin(pkin(5));
	t258 = sin(qJ(2));
	t287 = t254 * t258;
	t259 = sin(qJ(1));
	t286 = t254 * t259;
	t285 = t254 * t260;
	t262 = cos(qJ(1));
	t284 = t254 * t262;
	t283 = t259 * t258;
	t261 = cos(qJ(2));
	t282 = t259 * t261;
	t281 = t262 * t258;
	t280 = t262 * t261;
	t279 = qJD(1) * t259;
	t278 = qJD(1) * t262;
	t277 = qJD(2) * t258;
	t276 = qJD(2) * t261;
	t255 = cos(pkin(5));
	t241 = t255 * t281 + t282;
	t275 = qJD(3) * t241;
	t274 = qJD(3) * t262;
	t273 = t255 * t283;
	t272 = t255 * t280;
	t271 = t254 * t279;
	t270 = t254 * t278;
	t269 = t254 * t274;
	t268 = qJD(2) * t255 + qJD(1);
	t250 = pkin(2) + t289;
	t267 = t252 * r_i_i_C(1) + t250 - t290;
	t266 = t255 * t282 + t281;
	t265 = t271 - t275;
	t244 = t252 * t269;
	t243 = -t273 + t280;
	t240 = -t272 + t283;
	t239 = -qJD(1) * t273 - t259 * t277 + t268 * t280;
	t238 = t266 * qJD(1) + t241 * qJD(2);
	t237 = t241 * qJD(1) + t266 * qJD(2);
	t236 = -qJD(1) * t272 - t262 * t276 + t268 * t283;
	t235 = t251 * t270 - t237 * t252 + (-t243 * t251 + t252 * t286) * qJD(3);
	t234 = t252 * t270 + t237 * t251 + (-t243 * t252 - t251 * t286) * qJD(3);
	t1 = [(-t239 * t252 + t251 * t275 + t244) * r_i_i_C(1) + (t239 * t251 + t252 * t275) * r_i_i_C(2) - t239 * t250 + t275 * t291 - t240 * qJD(4) - pkin(1) * t278 - t288 * t238 + ((t289 - t290) * t274 + (-pkin(7) - t264) * t279) * t254, t243 * qJD(4) + t267 * t236 - t288 * t237 + t263 * t266, t234 * r_i_i_C(1) - t235 * r_i_i_C(2) + (t260 * t270 + t237 * t257 + (-t243 * t260 - t257 * t286) * qJD(3)) * pkin(3), -t236, 0; t235 * r_i_i_C(1) + t234 * r_i_i_C(2) + t266 * qJD(4) - t237 * t250 - t288 * t236 + (-pkin(1) * t259 + pkin(7) * t284) * qJD(1) + (t257 * t270 + (-t243 * t257 + t259 * t285) * qJD(3)) * pkin(3), t241 * qJD(4) - t267 * t238 + t288 * t239 + t240 * t263, t244 * r_i_i_C(2) + (t265 * r_i_i_C(1) - t239 * r_i_i_C(2)) * t252 + ((-t239 + t269) * r_i_i_C(1) - t265 * r_i_i_C(2)) * t251 + (t260 * t271 - t239 * t257 + (-t241 * t260 + t257 * t284) * qJD(3)) * pkin(3), t238, 0; 0, (qJD(4) * t258 - t261 * t263 + (-t267 * t258 + t288 * t261) * qJD(2)) * t254, -t264 * t254 * t276 + ((-t251 * t255 - t252 * t287) * r_i_i_C(1) + (t251 * t287 - t252 * t255) * r_i_i_C(2) + (-t255 * t257 - t258 * t285) * pkin(3)) * qJD(3), t254 * t277, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:30
	% EndTime: 2019-12-31 21:31:31
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (686->123), mult. (1483->200), div. (0->0), fcn. (1476->12), ass. (0->71)
	t393 = sin(qJ(1));
	t388 = cos(pkin(5));
	t409 = qJD(2) * t388 + qJD(1);
	t392 = sin(qJ(2));
	t427 = t392 * t393;
	t416 = t388 * t427;
	t422 = qJD(2) * t392;
	t396 = cos(qJ(2));
	t397 = cos(qJ(1));
	t424 = t396 * t397;
	t362 = -qJD(1) * t416 - t393 * t422 + t409 * t424;
	t425 = t393 * t396;
	t426 = t392 * t397;
	t373 = t388 * t426 + t425;
	t386 = qJ(3) + pkin(10);
	t384 = sin(t386);
	t385 = cos(t386);
	t387 = sin(pkin(5));
	t423 = qJD(1) * t387;
	t414 = t393 * t423;
	t428 = t387 * t397;
	t417 = t385 * t428;
	t356 = (-qJD(3) * t373 + t414) * t384 - qJD(3) * t417 + t362 * t385;
	t374 = t388 * t425 + t426;
	t361 = t374 * qJD(1) + t373 * qJD(2);
	t390 = sin(qJ(5));
	t394 = cos(qJ(5));
	t447 = t356 * t390 - t361 * t394;
	t446 = -t356 * t394 - t361 * t390;
	t391 = sin(qJ(3));
	t407 = r_i_i_C(1) * t394 - r_i_i_C(2) * t390;
	t405 = pkin(4) + t407;
	t440 = pkin(9) + r_i_i_C(3);
	t445 = (t391 * pkin(3) + t405 * t384 - t440 * t385) * qJD(3);
	t367 = -t373 * t385 + t384 * t428;
	t415 = t388 * t424;
	t372 = -t415 + t427;
	t444 = -t367 * t390 - t372 * t394;
	t443 = t367 * t394 - t372 * t390;
	t406 = r_i_i_C(1) * t390 + r_i_i_C(2) * t394;
	t395 = cos(qJ(3));
	t383 = pkin(3) * t395 + pkin(2);
	t442 = t440 * t384 + t405 * t385 + t383;
	t432 = t387 * t392;
	t431 = t387 * t393;
	t430 = t387 * t395;
	t429 = t387 * t396;
	t421 = qJD(2) * t396;
	t420 = qJD(5) * t385;
	t419 = qJD(5) * t390;
	t418 = qJD(5) * t394;
	t413 = t397 * t423;
	t412 = t387 * t421;
	t411 = t387 * t422;
	t375 = -t416 + t424;
	t404 = -t375 * t384 + t385 * t431;
	t369 = t375 * t385 + t384 * t431;
	t371 = t384 * t388 + t385 * t432;
	t403 = -t384 * t432 + t385 * t388;
	t402 = qJD(5) * t406;
	t399 = t367 * qJD(3) - t362 * t384 + t385 * t414;
	t398 = t406 * t420 + t445;
	t389 = -qJ(4) - pkin(8);
	t364 = t403 * qJD(3) + t385 * t412;
	t360 = t373 * qJD(1) + t374 * qJD(2);
	t359 = -qJD(1) * t415 - t397 * t421 + t409 * t427;
	t354 = t404 * qJD(3) - t360 * t385 + t384 * t413;
	t353 = t369 * qJD(3) - t360 * t384 - t385 * t413;
	t352 = t354 * t394 - t359 * t390 + (-t369 * t390 + t374 * t394) * qJD(5);
	t351 = -t354 * t390 - t359 * t394 + (-t369 * t394 - t374 * t390) * qJD(5);
	t1 = [t446 * r_i_i_C(1) + t447 * r_i_i_C(2) - t356 * pkin(4) - t362 * t383 + t361 * t389 - t372 * qJD(4) + t440 * t399 + (t444 * r_i_i_C(1) - t443 * r_i_i_C(2)) * qJD(5) + (-pkin(1) * t397 - pkin(7) * t431) * qJD(1) + (-t391 * t414 + (t373 * t391 + t395 * t428) * qJD(3)) * pkin(3), (-t360 * t390 + t375 * t418) * r_i_i_C(1) + (-t360 * t394 - t375 * t419) * r_i_i_C(2) + t360 * t389 + t375 * qJD(4) + t442 * t359 + t398 * t374, t440 * t354 - t404 * t402 - t405 * t353 + (t395 * t413 + t360 * t391 + (-t375 * t395 - t391 * t431) * qJD(3)) * pkin(3), -t359, r_i_i_C(1) * t351 - r_i_i_C(2) * t352; t354 * pkin(4) + t352 * r_i_i_C(1) + t351 * r_i_i_C(2) + t374 * qJD(4) + t359 * t389 - t360 * t383 + t440 * t353 + (-pkin(1) * t393 + pkin(7) * t428) * qJD(1) + (t391 * t413 + (-t375 * t391 + t393 * t430) * qJD(3)) * pkin(3), (t362 * t390 + t373 * t418) * r_i_i_C(1) + (t362 * t394 - t373 * t419) * r_i_i_C(2) - t362 * t389 + t373 * qJD(4) - t442 * t361 + t398 * t372, t440 * t356 - (-t373 * t384 - t417) * t402 + t405 * t399 + (t395 * t414 - t362 * t391 + (-t373 * t395 + t391 * t428) * qJD(3)) * pkin(3), t361, -t447 * r_i_i_C(1) + t446 * r_i_i_C(2) + (t443 * r_i_i_C(1) + t444 * r_i_i_C(2)) * qJD(5); 0, ((-qJD(2) * t442 + t407 * qJD(5) + qJD(4)) * t392 + (-qJD(2) * t389 - t445 + t406 * (qJD(2) - t420)) * t396) * t387, t440 * t364 - t403 * t402 + t405 * (-t371 * qJD(3) - t384 * t412) + (-t391 * t412 + (-t388 * t391 - t392 * t430) * qJD(3)) * pkin(3), t411, (-t364 * t390 + t394 * t411) * r_i_i_C(1) + (-t364 * t394 - t390 * t411) * r_i_i_C(2) + ((-t371 * t394 + t390 * t429) * r_i_i_C(1) + (t371 * t390 + t394 * t429) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end