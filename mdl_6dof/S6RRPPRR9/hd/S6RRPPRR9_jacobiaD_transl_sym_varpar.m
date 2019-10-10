% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:18
	% EndTime: 2019-10-10 09:50:18
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (84->29), mult. (244->45), div. (0->0), fcn. (222->6), ass. (0->24)
	t175 = sin(pkin(6));
	t194 = t175 * (pkin(8) + r_i_i_C(1));
	t192 = r_i_i_C(2) - pkin(2);
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
	t1 = [-(-t184 + t190) * qJD(3) + t192 * t170 - t191 * t169 + (-t180 * pkin(1) - t178 * t194) * qJD(1), -(t185 - t187) * qJD(3) - t191 * t168 - t192 * t167, -t167, 0, 0, 0; t182 * qJD(3) + t192 * t168 - t191 * t167 + (-t178 * pkin(1) + t180 * t194) * qJD(1), t181 * qJD(3) + t192 * t169 + t191 * t170, t169, 0, 0, 0; 0, (t177 * qJD(3) + (t192 * t177 + t191 * t179) * qJD(2)) * t175, t175 * t186, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:18
	% EndTime: 2019-10-10 09:50:18
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (118->36), mult. (337->52), div. (0->0), fcn. (312->6), ass. (0->27)
	t175 = sin(pkin(6));
	t193 = t175 * (pkin(3) + pkin(8) + r_i_i_C(1));
	t192 = r_i_i_C(2) + qJ(3);
	t177 = sin(qJ(2));
	t178 = sin(qJ(1));
	t191 = t178 * t177;
	t179 = cos(qJ(2));
	t190 = t178 * t179;
	t180 = cos(qJ(1));
	t189 = t180 * t177;
	t188 = t180 * t179;
	t187 = qJD(2) * t177;
	t186 = qJD(2) * t179;
	t184 = pkin(2) + r_i_i_C(3) + qJ(4);
	t176 = cos(pkin(6));
	t183 = t176 * t191;
	t182 = t176 * t188;
	t181 = qJD(2) * t176 + qJD(1);
	t169 = t176 * t190 + t189;
	t168 = t176 * t189 + t190;
	t170 = t183 - t188;
	t167 = -t182 + t191;
	t166 = -qJD(1) * t183 - t178 * t187 + t181 * t188;
	t165 = t169 * qJD(1) + t168 * qJD(2);
	t164 = t168 * qJD(1) + t169 * qJD(2);
	t163 = -qJD(1) * t182 - t180 * t186 + t181 * t191;
	t1 = [-t167 * qJD(3) - t168 * qJD(4) - t192 * t165 - t184 * t166 + (-pkin(1) * t180 - t178 * t193) * qJD(1), -t170 * qJD(3) - t169 * qJD(4) + t184 * t163 - t192 * t164, -t163, -t164, 0, 0; t169 * qJD(3) - t170 * qJD(4) - t192 * t163 - t184 * t164 + (-pkin(1) * t178 + t180 * t193) * qJD(1), t168 * qJD(3) - t167 * qJD(4) - t184 * t165 + t192 * t166, t165, t166, 0, 0; 0, (qJD(3) * t177 + qJD(4) * t179 + (-t184 * t177 + t192 * t179) * qJD(2)) * t175, t175 * t187, t175 * t186, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:18
	% EndTime: 2019-10-10 09:50:18
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (213->61), mult. (627->100), div. (0->0), fcn. (592->8), ass. (0->43)
	t248 = sin(qJ(5));
	t251 = cos(qJ(5));
	t257 = t251 * r_i_i_C(1) - t248 * r_i_i_C(2);
	t278 = t257 * qJD(5) + qJD(4);
	t277 = -pkin(2) - qJ(4);
	t246 = sin(pkin(6));
	t276 = t246 * t248;
	t275 = t246 * t251;
	t253 = cos(qJ(1));
	t274 = t246 * t253;
	t249 = sin(qJ(2));
	t250 = sin(qJ(1));
	t273 = t249 * t250;
	t272 = t249 * t253;
	t252 = cos(qJ(2));
	t271 = t250 * t252;
	t270 = t252 * t253;
	t269 = qJD(1) * t250;
	t268 = qJD(1) * t253;
	t267 = qJD(2) * t249;
	t266 = qJD(2) * t252;
	t265 = pkin(3) + pkin(4) + pkin(8);
	t264 = r_i_i_C(3) + pkin(9) - qJ(3);
	t247 = cos(pkin(6));
	t263 = t247 * t273;
	t262 = t247 * t270;
	t261 = t246 * t269;
	t260 = t246 * t268;
	t259 = t246 * t266;
	t258 = qJD(2) * t247 + qJD(1);
	t256 = -t248 * r_i_i_C(1) - t251 * r_i_i_C(2);
	t238 = t247 * t271 + t272;
	t237 = t247 * t272 + t271;
	t255 = -t256 - t277;
	t239 = -t263 + t270;
	t236 = -t262 + t273;
	t235 = -qJD(1) * t263 - t250 * t267 + t258 * t270;
	t234 = t238 * qJD(1) + t237 * qJD(2);
	t233 = t237 * qJD(1) + t238 * qJD(2);
	t232 = -qJD(1) * t262 - t253 * t266 + t258 * t273;
	t231 = t251 * t260 - t233 * t248 + (t239 * t251 - t250 * t276) * qJD(5);
	t230 = -t248 * t260 - t233 * t251 + (-t239 * t248 - t250 * t275) * qJD(5);
	t1 = [-pkin(1) * t268 - t236 * qJD(3) - t278 * t237 - t255 * t235 + t264 * t234 + (t256 * t253 * qJD(5) + (-t257 - t265) * t269) * t246, qJD(3) * t239 + t255 * t232 + t264 * t233 - t238 * t278, -t232, -t233, r_i_i_C(1) * t230 - r_i_i_C(2) * t231, 0; t231 * r_i_i_C(1) + t230 * r_i_i_C(2) + t238 * qJD(3) + t239 * qJD(4) + t277 * t233 + t264 * t232 + (-pkin(1) * t250 + t265 * t274) * qJD(1), qJD(3) * t237 - t255 * t234 - t264 * t235 - t236 * t278, t234, t235, (t235 * t251 - t248 * t261) * r_i_i_C(1) + (-t235 * t248 - t251 * t261) * r_i_i_C(2) + ((-t237 * t248 + t251 * t274) * r_i_i_C(1) + (-t237 * t251 - t248 * t274) * r_i_i_C(2)) * qJD(5), 0; 0, (qJD(3) * t249 + t278 * t252 + (-t255 * t249 - t264 * t252) * qJD(2)) * t246, t246 * t267, t259, t257 * t259 + ((-t247 * t251 - t249 * t276) * r_i_i_C(1) + (t247 * t248 - t249 * t275) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:20
	% EndTime: 2019-10-10 09:50:21
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (512->99), mult. (1517->159), div. (0->0), fcn. (1514->10), ass. (0->65)
	t390 = sin(qJ(2));
	t391 = sin(qJ(1));
	t394 = cos(qJ(2));
	t395 = cos(qJ(1));
	t433 = cos(pkin(6));
	t414 = t395 * t433;
	t375 = t390 * t414 + t391 * t394;
	t389 = sin(qJ(5));
	t393 = cos(qJ(5));
	t387 = sin(pkin(6));
	t427 = t387 * t395;
	t370 = -t375 * t389 + t393 * t427;
	t410 = t394 * t414;
	t426 = t390 * t391;
	t374 = -t410 + t426;
	t388 = sin(qJ(6));
	t392 = cos(qJ(6));
	t442 = -t370 * t388 + t374 * t392;
	t441 = t370 * t392 + t374 * t388;
	t408 = t388 * r_i_i_C(1) + t392 * r_i_i_C(2);
	t403 = qJD(6) * t408;
	t409 = -t392 * r_i_i_C(1) + t388 * r_i_i_C(2);
	t406 = pkin(5) - t409;
	t436 = pkin(10) + r_i_i_C(3);
	t396 = -(t436 * t389 + t406 * t393) * qJD(5) - qJD(4) + t389 * t403;
	t407 = qJD(2) * t433 + qJD(1);
	t415 = t391 * t433;
	t411 = t390 * t415;
	t423 = qJD(2) * t390;
	t425 = t394 * t395;
	t364 = -qJD(1) * t411 - t391 * t423 + t407 * t425;
	t429 = t387 * t393;
	t418 = qJD(5) * t429;
	t424 = qJD(1) * t387;
	t420 = t391 * t424;
	t439 = (qJD(5) * t375 + t420) * t389 - t364 * t393 - t395 * t418;
	t405 = t375 * t393 + t389 * t427;
	t356 = t405 * qJD(5) + t364 * t389 + t393 * t420;
	t435 = -pkin(2) - qJ(4);
	t437 = t406 * t389 - t436 * t393 - t435;
	t434 = pkin(9) - qJ(3);
	t430 = t387 * t391;
	t428 = t387 * t394;
	t422 = qJD(2) * t394;
	t421 = pkin(3) + pkin(4) + pkin(8);
	t419 = t395 * t424;
	t417 = t387 * t423;
	t416 = t387 * t422;
	t413 = t433 * t389;
	t377 = -t411 + t425;
	t368 = t377 * t389 + t391 * t429;
	t404 = t377 * t393 - t389 * t430;
	t402 = t408 + t434;
	t376 = t395 * t390 + t394 * t415;
	t399 = -t387 * t390 * t389 - t433 * t393;
	t398 = t409 * qJD(6) + qJD(3);
	t365 = qJD(5) * t413 - t389 * t416 - t390 * t418;
	t363 = t376 * qJD(1) + t375 * qJD(2);
	t362 = t375 * qJD(1) + t376 * qJD(2);
	t361 = -qJD(1) * t410 - t395 * t422 + t407 * t426;
	t359 = t404 * qJD(5) - t362 * t389 + t393 * t419;
	t358 = t368 * qJD(5) + t362 * t393 + t389 * t419;
	t353 = t359 * t392 + t361 * t388 + (-t368 * t388 - t376 * t392) * qJD(6);
	t352 = -t359 * t388 + t361 * t392 + (-t368 * t392 + t376 * t388) * qJD(6);
	t1 = [-t374 * qJD(3) - t375 * qJD(4) + t435 * t364 - t406 * t356 + t402 * t363 - t436 * t439 + (t442 * r_i_i_C(1) - t441 * r_i_i_C(2)) * qJD(6) + (-t395 * pkin(1) - t421 * t430) * qJD(1), t437 * t361 + t402 * t362 + t396 * t376 + t398 * t377, -t361, -t362, -t406 * t358 + t436 * t359 - t404 * t403, r_i_i_C(1) * t352 - r_i_i_C(2) * t353; t359 * pkin(5) + t353 * r_i_i_C(1) + t352 * r_i_i_C(2) + t376 * qJD(3) + t377 * qJD(4) + t435 * t362 + t434 * t361 + t436 * t358 + (-pkin(1) * t391 + t421 * t427) * qJD(1), -t363 * t437 - t402 * t364 + t396 * t374 + t398 * t375, t363, t364, t436 * t356 - t405 * t403 - t406 * t439, (-t356 * t388 - t363 * t392) * r_i_i_C(1) + (-t356 * t392 + t363 * t388) * r_i_i_C(2) + (t441 * r_i_i_C(1) + t442 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t437 + t398) * t390 + (-t402 * qJD(2) - t396) * t394) * t387, t417, t416, -t436 * t365 - (t390 * t429 - t413) * t403 + t406 * (t399 * qJD(5) + t393 * t416), (t365 * t388 - t392 * t417) * r_i_i_C(1) + (t365 * t392 + t388 * t417) * r_i_i_C(2) + ((-t388 * t428 + t392 * t399) * r_i_i_C(1) + (-t388 * t399 - t392 * t428) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end