% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(10));
	t50 = sin(pkin(10));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (22->14), mult. (83->33), div. (0->0), fcn. (72->8), ass. (0->16)
	t89 = sin(pkin(11));
	t92 = cos(pkin(11));
	t95 = sin(qJ(2));
	t96 = cos(qJ(2));
	t98 = t89 * t96 + t92 * t95;
	t88 = t98 * qJD(2);
	t94 = cos(pkin(6));
	t100 = t94 * t95;
	t99 = pkin(2) * qJD(2);
	t97 = t89 * t95 - t92 * t96;
	t87 = t97 * qJD(2);
	t93 = cos(pkin(10));
	t90 = sin(pkin(10));
	t86 = t94 * t88;
	t85 = t94 * t87;
	t1 = [0, (t90 * t86 + t93 * t87) * r_i_i_C(1) + (-t90 * t85 + t93 * t88) * r_i_i_C(2) + (t90 * t100 - t93 * t96) * t99, 0, 0, 0, 0; 0, (-t93 * t86 + t90 * t87) * r_i_i_C(1) + (t93 * t85 + t90 * t88) * r_i_i_C(2) + (-t93 * t100 - t90 * t96) * t99, 0, 0, 0, 0; 0, (-t95 * pkin(2) - t98 * r_i_i_C(1) + t97 * r_i_i_C(2)) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:30:00
	% EndTime: 2019-10-09 21:30:00
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (115->39), mult. (386->84), div. (0->0), fcn. (382->10), ass. (0->37)
	t241 = sin(pkin(11));
	t244 = cos(pkin(11));
	t250 = cos(qJ(2));
	t261 = qJD(2) * t250;
	t248 = sin(qJ(2));
	t262 = qJD(2) * t248;
	t268 = t241 * t262 - t244 * t261;
	t267 = -pkin(8) - r_i_i_C(3);
	t266 = pkin(2) * qJD(2);
	t243 = sin(pkin(6));
	t247 = sin(qJ(4));
	t265 = t243 * t247;
	t249 = cos(qJ(4));
	t264 = t243 * t249;
	t246 = cos(pkin(6));
	t263 = t246 * t248;
	t258 = r_i_i_C(1) * t247 + r_i_i_C(2) * t249;
	t227 = t268 * t246;
	t234 = -t241 * t261 - t244 * t262;
	t242 = sin(pkin(10));
	t245 = cos(pkin(10));
	t257 = t227 * t245 - t234 * t242;
	t256 = t227 * t242 + t234 * t245;
	t255 = t241 * t250 + t248 * t244;
	t254 = t248 * t241 - t244 * t250;
	t253 = r_i_i_C(1) * t249 - r_i_i_C(2) * t247 + pkin(3);
	t252 = qJD(4) * t258;
	t251 = qJD(2) * t255;
	t233 = t254 * qJD(2);
	t232 = t255 * t246;
	t231 = t254 * t246;
	t230 = t255 * t243;
	t228 = t246 * t251;
	t225 = t268 * t243;
	t224 = -t232 * t242 - t245 * t254;
	t222 = t232 * t245 - t242 * t254;
	t1 = [0, -t267 * t256 - (t231 * t242 - t245 * t255) * t252 + (t242 * t263 - t245 * t250) * t266 + t253 * (t228 * t242 + t233 * t245), 0, -t258 * t256 + ((-t224 * t249 - t242 * t265) * r_i_i_C(1) + (t224 * t247 - t242 * t264) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t267 * t257 - (-t231 * t245 - t242 * t255) * t252 + (-t242 * t250 - t245 * t263) * t266 + t253 * (-t228 * t245 + t233 * t242), 0, t258 * t257 + ((-t222 * t249 + t245 * t265) * r_i_i_C(1) + (t222 * t247 + t245 * t264) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t267 * t225 + (-pkin(2) * t262 - t251 * t253 + t252 * t254) * t243, 0, t258 * t225 + ((-t230 * t249 - t246 * t247) * r_i_i_C(1) + (t230 * t247 - t246 * t249) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:30:00
	% EndTime: 2019-10-09 21:30:00
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (201->53), mult. (547->99), div. (0->0), fcn. (552->12), ass. (0->43)
	t261 = cos(pkin(6));
	t256 = sin(pkin(11));
	t259 = cos(pkin(11));
	t264 = sin(qJ(2));
	t266 = cos(qJ(2));
	t270 = t266 * t256 + t264 * t259;
	t242 = t270 * t261;
	t276 = qJD(2) * t266;
	t277 = qJD(2) * t264;
	t284 = t256 * t277 - t259 * t276;
	t283 = -r_i_i_C(3) - qJ(5) - pkin(8);
	t282 = pkin(2) * qJD(2);
	t257 = sin(pkin(10));
	t258 = sin(pkin(6));
	t281 = t257 * t258;
	t260 = cos(pkin(10));
	t280 = t258 * t260;
	t263 = sin(qJ(4));
	t279 = t258 * t263;
	t278 = t261 * t264;
	t237 = t284 * t261;
	t244 = -t256 * t276 - t259 * t277;
	t273 = t260 * t237 - t257 * t244;
	t272 = t257 * t237 + t260 * t244;
	t245 = t264 * t256 - t266 * t259;
	t232 = t260 * t242 - t257 * t245;
	t271 = t257 * t242 + t260 * t245;
	t255 = qJ(4) + pkin(12);
	t253 = sin(t255);
	t254 = cos(t255);
	t265 = cos(qJ(4));
	t269 = t265 * pkin(4) + t254 * r_i_i_C(1) - t253 * r_i_i_C(2) + pkin(3);
	t240 = t270 * t258;
	t268 = t263 * pkin(4) + t253 * r_i_i_C(1) + t254 * r_i_i_C(2);
	t267 = qJD(4) * t268;
	t243 = t245 * qJD(2);
	t241 = t245 * t261;
	t238 = qJD(2) * t242;
	t236 = qJD(2) * t240;
	t235 = t284 * t258;
	t228 = t257 * t238 + t260 * t243;
	t225 = -t260 * t238 + t257 * t243;
	t1 = [0, -t271 * qJD(5) - t283 * t272 + (t257 * t278 - t260 * t266) * t282 + t269 * t228 - (t257 * t241 - t260 * t270) * t267, 0, -t268 * t272 + ((-t253 * t281 + t254 * t271) * r_i_i_C(1) + (-t253 * t271 - t254 * t281) * r_i_i_C(2) + (-t257 * t279 + t265 * t271) * pkin(4)) * qJD(4), -t228, 0; 0, t232 * qJD(5) + t283 * t273 + (-t257 * t266 - t260 * t278) * t282 + t269 * t225 - (-t260 * t241 - t257 * t270) * t267, 0, t268 * t273 + ((-t232 * t254 + t253 * t280) * r_i_i_C(1) + (t232 * t253 + t254 * t280) * r_i_i_C(2) + (-t232 * t265 + t260 * t279) * pkin(4)) * qJD(4), -t225, 0; 0, t240 * qJD(5) + t283 * t235 - t269 * t236 + (-pkin(2) * t277 + t245 * t267) * t258, 0, t268 * t235 + ((-t240 * t254 - t253 * t261) * r_i_i_C(1) + (t240 * t253 - t254 * t261) * r_i_i_C(2) + (-t240 * t265 - t261 * t263) * pkin(4)) * qJD(4), t236, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:30:01
	% EndTime: 2019-10-09 21:30:02
	% DurationCPUTime: 0.83s
	% Computational Cost: add. (608->99), mult. (1513->178), div. (0->0), fcn. (1626->14), ass. (0->61)
	t390 = qJ(4) + pkin(12);
	t388 = sin(t390);
	t389 = cos(t390);
	t398 = sin(qJ(4));
	t397 = sin(qJ(6));
	t400 = cos(qJ(6));
	t408 = (t397 * r_i_i_C(1) + t400 * r_i_i_C(2)) * qJD(6);
	t413 = t400 * r_i_i_C(1) - t397 * r_i_i_C(2) + pkin(5);
	t432 = pkin(9) + r_i_i_C(3);
	t436 = (t398 * pkin(4) + t413 * t388 - t432 * t389) * qJD(4) + t389 * t408;
	t395 = cos(pkin(6));
	t391 = sin(pkin(11));
	t399 = sin(qJ(2));
	t428 = cos(pkin(11));
	t431 = cos(qJ(2));
	t407 = t431 * t391 + t399 * t428;
	t376 = t407 * t395;
	t416 = t431 * t428;
	t422 = qJD(2) * t399;
	t434 = -qJD(2) * t416 + t391 * t422;
	t406 = -t399 * t391 + t416;
	t401 = cos(qJ(4));
	t403 = t401 * pkin(4) + t432 * t388 + t413 * t389 + pkin(3);
	t429 = pkin(2) * qJD(2);
	t392 = sin(pkin(10));
	t393 = sin(pkin(6));
	t427 = t392 * t393;
	t394 = cos(pkin(10));
	t426 = t393 * t394;
	t425 = t393 * t398;
	t424 = t395 * t399;
	t421 = qJD(6) * t397;
	t420 = qJD(6) * t400;
	t373 = t434 * t395;
	t378 = t407 * qJD(2);
	t355 = t394 * t373 + t392 * t378;
	t357 = t392 * t373 - t394 * t378;
	t375 = t407 * t393;
	t366 = t375 * t389 + t395 * t388;
	t414 = -t375 * t388 + t395 * t389;
	t361 = t394 * t376 + t392 * t406;
	t362 = t392 * t376 - t394 * t406;
	t411 = -t361 * t388 - t389 * t426;
	t410 = -t361 * t389 + t388 * t426;
	t409 = t362 * t388 + t389 * t427;
	t352 = -t362 * t389 + t388 * t427;
	t405 = t406 * t395;
	t404 = qJD(2) * t376;
	t396 = -qJ(5) - pkin(8);
	t377 = t406 * qJD(2);
	t374 = t406 * t393;
	t372 = qJD(2) * t375;
	t371 = t434 * t393;
	t363 = -t392 * t405 - t394 * t407;
	t360 = -t392 * t407 + t394 * t405;
	t356 = -t394 * t377 + t392 * t404;
	t353 = -t392 * t377 - t394 * t404;
	t348 = t414 * qJD(4) - t371 * t389;
	t346 = t409 * qJD(4) + t357 * t389;
	t344 = t411 * qJD(4) - t355 * t389;
	t1 = [0, (t357 * t397 - t362 * t420) * r_i_i_C(1) + (t357 * t400 + t362 * t421) * r_i_i_C(2) - t357 * t396 - t362 * qJD(5) + (t392 * t424 - t431 * t394) * t429 + t403 * t356 - t436 * t363, 0, t432 * t346 - t409 * t408 + t413 * (-t352 * qJD(4) - t357 * t388) + (-t357 * t398 + (t362 * t401 - t392 * t425) * qJD(4)) * pkin(4), -t356, (-t346 * t397 - t356 * t400) * r_i_i_C(1) + (-t346 * t400 + t356 * t397) * r_i_i_C(2) + ((-t352 * t400 + t363 * t397) * r_i_i_C(1) + (t352 * t397 + t363 * t400) * r_i_i_C(2)) * qJD(6); 0, (-t355 * t397 + t361 * t420) * r_i_i_C(1) + (-t355 * t400 - t361 * t421) * r_i_i_C(2) + t355 * t396 + t361 * qJD(5) + (-t431 * t392 - t394 * t424) * t429 + t403 * t353 - t436 * t360, 0, t432 * t344 - t411 * t408 + t413 * (t410 * qJD(4) + t355 * t388) + (t355 * t398 + (-t361 * t401 + t394 * t425) * qJD(4)) * pkin(4), -t353, (-t344 * t397 - t353 * t400) * r_i_i_C(1) + (-t344 * t400 + t353 * t397) * r_i_i_C(2) + ((t360 * t397 + t400 * t410) * r_i_i_C(1) + (t360 * t400 - t397 * t410) * r_i_i_C(2)) * qJD(6); 0, (-t371 * t397 + t375 * t420) * r_i_i_C(1) + (-t371 * t400 - t375 * t421) * r_i_i_C(2) + t371 * t396 + t375 * qJD(5) - t393 * pkin(2) * t422 - t403 * t372 - t436 * t374, 0, t432 * t348 - t414 * t408 + t413 * (-t366 * qJD(4) + t371 * t388) + (t371 * t398 + (-t375 * t401 - t395 * t398) * qJD(4)) * pkin(4), t372, (-t348 * t397 + t372 * t400) * r_i_i_C(1) + (-t348 * t400 - t372 * t397) * r_i_i_C(2) + ((-t366 * t400 + t374 * t397) * r_i_i_C(1) + (t366 * t397 + t374 * t400) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end