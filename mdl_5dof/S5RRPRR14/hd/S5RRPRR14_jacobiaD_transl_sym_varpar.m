% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR14
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR14_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR14_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR14_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:24:51
	% EndTime: 2019-12-29 19:24:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:24:59
	% EndTime: 2019-12-29 19:24:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:24:59
	% EndTime: 2019-12-29 19:24:59
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-12-29 19:24:55
	% EndTime: 2019-12-29 19:24:55
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (103->31), mult. (307->49), div. (0->0), fcn. (282->8), ass. (0->26)
	t199 = sin(pkin(10));
	t200 = sin(pkin(5));
	t201 = cos(pkin(10));
	t220 = t200 * (r_i_i_C(1) * t199 + r_i_i_C(2) * t201 + pkin(7));
	t219 = r_i_i_C(3) + qJ(3);
	t203 = sin(qJ(2));
	t204 = sin(qJ(1));
	t218 = t204 * t203;
	t205 = cos(qJ(2));
	t217 = t204 * t205;
	t206 = cos(qJ(1));
	t216 = t206 * t203;
	t215 = t206 * t205;
	t214 = qJD(2) * t203;
	t202 = cos(pkin(5));
	t213 = t202 * t218;
	t212 = t202 * t215;
	t211 = qJD(2) * t202 + qJD(1);
	t210 = t201 * r_i_i_C(1) - t199 * r_i_i_C(2) + pkin(2);
	t208 = t202 * t217 + t216;
	t207 = t202 * t216 + t217;
	t194 = -qJD(1) * t213 - t204 * t214 + t211 * t215;
	t193 = t208 * qJD(1) + t207 * qJD(2);
	t192 = t207 * qJD(1) + t208 * qJD(2);
	t191 = -qJD(1) * t212 - qJD(2) * t215 + t211 * t218;
	t1 = [-(-t212 + t218) * qJD(3) - t219 * t193 - t210 * t194 + (-t206 * pkin(1) - t204 * t220) * qJD(1), -(t213 - t215) * qJD(3) - t219 * t192 + t210 * t191, -t191, 0, 0; t208 * qJD(3) - t219 * t191 - t210 * t192 + (-t204 * pkin(1) + t206 * t220) * qJD(1), t207 * qJD(3) - t210 * t193 + t219 * t194, t193, 0, 0; 0, (t203 * qJD(3) + (-t210 * t203 + t219 * t205) * qJD(2)) * t200, t200 * t214, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:24:53
	% EndTime: 2019-12-29 19:24:54
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (222->61), mult. (494->100), div. (0->0), fcn. (467->10), ass. (0->45)
	t258 = sin(qJ(1));
	t255 = cos(pkin(5));
	t266 = qJD(2) * t255 + qJD(1);
	t257 = sin(qJ(2));
	t279 = t258 * t257;
	t271 = t255 * t279;
	t274 = qJD(2) * t257;
	t259 = cos(qJ(2));
	t260 = cos(qJ(1));
	t276 = t260 * t259;
	t238 = -qJD(1) * t271 - t258 * t274 + t266 * t276;
	t254 = sin(pkin(5));
	t280 = t254 * t260;
	t267 = qJD(4) * t280;
	t284 = t238 - t267;
	t252 = pkin(10) + qJ(4);
	t250 = sin(t252);
	t251 = cos(t252);
	t265 = r_i_i_C(1) * t250 + r_i_i_C(2) * t251;
	t262 = qJD(4) * t265;
	t283 = r_i_i_C(3) + pkin(8) + qJ(3);
	t282 = t254 * t257;
	t281 = t254 * t258;
	t278 = t258 * t259;
	t277 = t260 * t257;
	t275 = qJD(1) * t254;
	t273 = qJD(2) * t259;
	t240 = t255 * t277 + t278;
	t272 = qJD(4) * t240;
	t270 = t255 * t276;
	t269 = pkin(3) * sin(pkin(10)) + pkin(7);
	t268 = t260 * t275;
	t249 = cos(pkin(10)) * pkin(3) + pkin(2);
	t264 = t251 * r_i_i_C(1) - t250 * r_i_i_C(2) + t249;
	t263 = t255 * t278 + t277;
	t261 = t258 * t275 - t272;
	t243 = t251 * t267;
	t242 = -t271 + t276;
	t239 = -t270 + t279;
	t237 = t263 * qJD(1) + t240 * qJD(2);
	t236 = t240 * qJD(1) + t263 * qJD(2);
	t235 = -qJD(1) * t270 - t260 * t273 + t266 * t279;
	t234 = t250 * t268 - t236 * t251 + (-t242 * t250 + t251 * t281) * qJD(4);
	t233 = t251 * t268 + t236 * t250 + (-t242 * t251 - t250 * t281) * qJD(4);
	t1 = [(-t238 * t251 + t250 * t272 + t243) * r_i_i_C(1) + (t284 * t250 + t251 * t272) * r_i_i_C(2) - t238 * t249 - t239 * qJD(3) - t283 * t237 + (-t260 * pkin(1) + (-t265 - t269) * t281) * qJD(1), t242 * qJD(3) + t264 * t235 - t283 * t236 + t262 * t263, -t235, t233 * r_i_i_C(1) - t234 * r_i_i_C(2), 0; t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t263 * qJD(3) - t236 * t249 - t283 * t235 + (-pkin(1) * t258 + t269 * t280) * qJD(1), t240 * qJD(3) - t264 * t237 + t283 * t238 + t239 * t262, t237, t243 * r_i_i_C(2) + (t261 * r_i_i_C(1) - t238 * r_i_i_C(2)) * t251 + (-t284 * r_i_i_C(1) - t261 * r_i_i_C(2)) * t250, 0; 0, (qJD(3) * t257 - t259 * t262 + (-t264 * t257 + t283 * t259) * qJD(2)) * t254, t254 * t274, -t265 * t254 * t273 + ((-t250 * t255 - t251 * t282) * r_i_i_C(1) + (t250 * t282 - t251 * t255) * r_i_i_C(2)) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:25:03
	% EndTime: 2019-12-29 19:25:04
	% DurationCPUTime: 1.17s
	% Computational Cost: add. (656->106), mult. (1384->172), div. (0->0), fcn. (1389->12), ass. (0->69)
	t392 = sin(qJ(1));
	t388 = cos(pkin(5));
	t407 = qJD(2) * t388 + qJD(1);
	t391 = sin(qJ(2));
	t426 = t392 * t391;
	t416 = t388 * t426;
	t421 = qJD(2) * t391;
	t394 = cos(qJ(2));
	t395 = cos(qJ(1));
	t423 = t395 * t394;
	t361 = -qJD(1) * t416 - t392 * t421 + t407 * t423;
	t424 = t395 * t391;
	t425 = t392 * t394;
	t372 = t388 * t424 + t425;
	t385 = pkin(10) + qJ(4);
	t383 = sin(t385);
	t384 = cos(t385);
	t387 = sin(pkin(5));
	t422 = qJD(1) * t387;
	t412 = t392 * t422;
	t427 = t387 * t395;
	t415 = t384 * t427;
	t355 = (-qJD(4) * t372 + t412) * t383 - qJD(4) * t415 + t361 * t384;
	t373 = t388 * t425 + t424;
	t360 = t373 * qJD(1) + t372 * qJD(2);
	t390 = sin(qJ(5));
	t393 = cos(qJ(5));
	t444 = t355 * t390 - t360 * t393;
	t443 = -t355 * t393 - t360 * t390;
	t405 = t393 * r_i_i_C(1) - t390 * r_i_i_C(2);
	t403 = pkin(4) + t405;
	t437 = r_i_i_C(3) + pkin(9);
	t442 = (t403 * t383 - t437 * t384) * qJD(4);
	t366 = -t372 * t384 + t383 * t427;
	t414 = t388 * t423;
	t371 = -t414 + t426;
	t441 = -t366 * t390 - t371 * t393;
	t440 = t366 * t393 - t371 * t390;
	t404 = t390 * r_i_i_C(1) + t393 * r_i_i_C(2);
	t382 = cos(pkin(10)) * pkin(3) + pkin(2);
	t438 = t437 * t383 + t403 * t384 + t382;
	t430 = t387 * t391;
	t429 = t387 * t392;
	t428 = t387 * t394;
	t420 = qJD(2) * t394;
	t419 = qJD(5) * t384;
	t418 = qJD(5) * t390;
	t417 = qJD(5) * t393;
	t413 = pkin(3) * sin(pkin(10)) + pkin(7);
	t411 = t395 * t422;
	t410 = t387 * t420;
	t409 = t387 * t421;
	t374 = -t416 + t423;
	t402 = -t374 * t383 + t384 * t429;
	t368 = t374 * t384 + t383 * t429;
	t370 = t388 * t383 + t384 * t430;
	t401 = -t383 * t430 + t388 * t384;
	t400 = qJD(5) * t404;
	t397 = t366 * qJD(4) - t361 * t383 + t384 * t412;
	t396 = t404 * t419 + t442;
	t389 = -pkin(8) - qJ(3);
	t363 = t401 * qJD(4) + t384 * t410;
	t359 = t372 * qJD(1) + t373 * qJD(2);
	t358 = -qJD(1) * t414 - t395 * t420 + t407 * t426;
	t353 = t402 * qJD(4) - t359 * t384 + t383 * t411;
	t352 = t368 * qJD(4) - t359 * t383 - t384 * t411;
	t351 = t353 * t393 - t358 * t390 + (-t368 * t390 + t373 * t393) * qJD(5);
	t350 = -t353 * t390 - t358 * t393 + (-t368 * t393 - t373 * t390) * qJD(5);
	t1 = [t443 * r_i_i_C(1) + t444 * r_i_i_C(2) - t355 * pkin(4) - t361 * t382 + t360 * t389 - t371 * qJD(3) + t437 * t397 + (t441 * r_i_i_C(1) - t440 * r_i_i_C(2)) * qJD(5) + (-t395 * pkin(1) - t413 * t429) * qJD(1), (-t359 * t390 + t374 * t417) * r_i_i_C(1) + (-t359 * t393 - t374 * t418) * r_i_i_C(2) + t359 * t389 + t374 * qJD(3) + t438 * t358 + t396 * t373, -t358, -t403 * t352 + t437 * t353 - t402 * t400, t350 * r_i_i_C(1) - t351 * r_i_i_C(2); t353 * pkin(4) + t351 * r_i_i_C(1) + t350 * r_i_i_C(2) + t373 * qJD(3) + t358 * t389 - t359 * t382 + t437 * t352 + (-pkin(1) * t392 + t413 * t427) * qJD(1), (t361 * t390 + t372 * t417) * r_i_i_C(1) + (t361 * t393 - t372 * t418) * r_i_i_C(2) - t361 * t389 + t372 * qJD(3) - t438 * t360 + t396 * t371, t360, t437 * t355 - (-t372 * t383 - t415) * t400 + t403 * t397, -t444 * r_i_i_C(1) + t443 * r_i_i_C(2) + (t440 * r_i_i_C(1) + t441 * r_i_i_C(2)) * qJD(5); 0, ((-qJD(2) * t438 + t405 * qJD(5) + qJD(3)) * t391 + (-qJD(2) * t389 - t442 + t404 * (qJD(2) - t419)) * t394) * t387, t409, t437 * t363 - t401 * t400 + t403 * (-t370 * qJD(4) - t383 * t410), (-t363 * t390 + t393 * t409) * r_i_i_C(1) + (-t363 * t393 - t390 * t409) * r_i_i_C(2) + ((-t370 * t393 + t390 * t428) * r_i_i_C(1) + (t370 * t390 + t393 * t428) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end