% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR16
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:30
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR16_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR16_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR16_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:31
	% EndTime: 2019-12-29 19:30:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:30
	% EndTime: 2019-12-29 19:30:30
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
	% StartTime: 2019-12-29 19:30:36
	% EndTime: 2019-12-29 19:30:36
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-12-29 19:30:31
	% EndTime: 2019-12-29 19:30:31
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (84->29), mult. (244->45), div. (0->0), fcn. (222->6), ass. (0->24)
	t175 = sin(pkin(5));
	t194 = t175 * (pkin(7) + r_i_i_C(1));
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
	t176 = cos(pkin(5));
	t185 = t176 * t190;
	t184 = t176 * t187;
	t183 = qJD(2) * t176 + qJD(1);
	t182 = t176 * t189 + t188;
	t181 = t176 * t188 + t189;
	t170 = -qJD(1) * t185 - t178 * t186 + t183 * t187;
	t169 = t182 * qJD(1) + t181 * qJD(2);
	t168 = t181 * qJD(1) + t182 * qJD(2);
	t167 = -qJD(1) * t184 - qJD(2) * t187 + t183 * t190;
	t1 = [-(-t184 + t190) * qJD(3) - t193 * t170 - t191 * t169 + (-t180 * pkin(1) - t178 * t194) * qJD(1), -(t185 - t187) * qJD(3) - t191 * t168 + t193 * t167, -t167, 0, 0; t182 * qJD(3) - t193 * t168 - t191 * t167 + (-t178 * pkin(1) + t180 * t194) * qJD(1), t181 * qJD(3) - t193 * t169 + t191 * t170, t169, 0, 0; 0, (t177 * qJD(3) + (-t193 * t177 + t191 * t179) * qJD(2)) * t175, t175 * t186, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:31
	% EndTime: 2019-12-29 19:30:32
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (179->58), mult. (534->98), div. (0->0), fcn. (502->8), ass. (0->43)
	t274 = pkin(3) + pkin(7);
	t242 = cos(pkin(5));
	t244 = sin(qJ(2));
	t248 = cos(qJ(1));
	t266 = t248 * t244;
	t245 = sin(qJ(1));
	t247 = cos(qJ(2));
	t267 = t245 * t247;
	t234 = t242 * t267 + t266;
	t250 = t242 * t266 + t267;
	t230 = t234 * qJD(1) + t250 * qJD(2);
	t243 = sin(qJ(4));
	t273 = t230 * t243;
	t246 = cos(qJ(4));
	t272 = t230 * t246;
	t241 = sin(pkin(5));
	t271 = t241 * t245;
	t270 = t241 * t247;
	t269 = t241 * t248;
	t268 = t245 * t244;
	t265 = t248 * t247;
	t264 = qJD(1) * t245;
	t263 = qJD(1) * t248;
	t262 = qJD(2) * t244;
	t258 = t242 * t265;
	t232 = -t258 + t268;
	t261 = qJD(4) * t232;
	t260 = -r_i_i_C(3) - pkin(8) - pkin(2);
	t259 = t242 * t268;
	t257 = t241 * t264;
	t256 = t241 * t263;
	t255 = t241 * t262;
	t254 = qJD(2) * t242 + qJD(1);
	t253 = r_i_i_C(1) * t246 - r_i_i_C(2) * t243;
	t252 = -t243 * r_i_i_C(1) - t246 * r_i_i_C(2);
	t251 = qJ(3) - t252;
	t249 = t253 * qJD(4) + qJD(3);
	t231 = -qJD(1) * t259 - t245 * t262 + t254 * t265;
	t229 = t250 * qJD(1) + t234 * qJD(2);
	t228 = -qJD(1) * t258 - qJD(2) * t265 + t254 * t268;
	t227 = t246 * t256 - t228 * t243 + (t234 * t246 - t243 * t271) * qJD(4);
	t226 = -t243 * t256 - t228 * t246 + (-t234 * t243 - t246 * t271) * qJD(4);
	t1 = [(-t246 * t261 - t273) * r_i_i_C(1) + (t243 * t261 - t272) * r_i_i_C(2) - t230 * qJ(3) - t232 * qJD(3) - pkin(1) * t263 + t260 * t231 + (t252 * t248 * qJD(4) + (-t253 - t274) * t264) * t241, t249 * (-t259 + t265) - t251 * t229 - t260 * t228, -t228, t226 * r_i_i_C(1) - t227 * r_i_i_C(2), 0; t227 * r_i_i_C(1) + t226 * r_i_i_C(2) - t228 * qJ(3) + t234 * qJD(3) + t260 * t229 + (-pkin(1) * t245 + t274 * t269) * qJD(1), t230 * t260 + t231 * t251 + t249 * t250, t230, (-t243 * t257 + t272) * r_i_i_C(1) + (-t246 * t257 - t273) * r_i_i_C(2) + ((-t232 * t243 + t246 * t269) * r_i_i_C(1) + (-t232 * t246 - t243 * t269) * r_i_i_C(2)) * qJD(4), 0; 0, (t249 * t244 + (t244 * t260 + t247 * t251) * qJD(2)) * t241, t255, t253 * t255 + ((-t242 * t246 + t243 * t270) * r_i_i_C(1) + (t242 * t243 + t246 * t270) * r_i_i_C(2)) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:34
	% EndTime: 2019-12-29 19:30:35
	% DurationCPUTime: 0.99s
	% Computational Cost: add. (478->92), mult. (1424->150), div. (0->0), fcn. (1424->10), ass. (0->61)
	t382 = cos(pkin(5));
	t389 = cos(qJ(2));
	t390 = cos(qJ(1));
	t415 = t390 * t389;
	t385 = sin(qJ(2));
	t386 = sin(qJ(1));
	t418 = t386 * t385;
	t372 = -t382 * t415 + t418;
	t384 = sin(qJ(4));
	t388 = cos(qJ(4));
	t381 = sin(pkin(5));
	t419 = t381 * t390;
	t411 = t388 * t419;
	t368 = -t372 * t384 + t411;
	t416 = t390 * t385;
	t417 = t386 * t389;
	t373 = t382 * t416 + t417;
	t383 = sin(qJ(5));
	t387 = cos(qJ(5));
	t430 = -t368 * t383 - t373 * t387;
	t429 = t368 * t387 - t373 * t383;
	t405 = t387 * r_i_i_C(1) - t383 * r_i_i_C(2);
	t396 = t405 * qJD(5);
	t374 = t382 * t417 + t416;
	t361 = t374 * qJD(1) + t373 * qJD(2);
	t414 = qJD(1) * t381;
	t410 = t386 * t414;
	t428 = (qJD(4) * t372 + t410) * t384 - qJD(4) * t411 - t361 * t388;
	t401 = t372 * t388 + t384 * t419;
	t354 = t401 * qJD(4) + t361 * t384 + t388 * t410;
	t403 = pkin(4) + t405;
	t425 = -r_i_i_C(3) - pkin(9);
	t392 = t403 * t384 + t425 * t388 + qJ(3);
	t427 = -pkin(2) - pkin(8);
	t426 = pkin(3) + pkin(7);
	t422 = t381 * t385;
	t421 = t381 * t386;
	t420 = t381 * t389;
	t413 = qJD(2) * t385;
	t412 = t382 * t418;
	t409 = t390 * t414;
	t408 = qJD(2) * t420;
	t407 = t381 * t413;
	t404 = -t383 * r_i_i_C(1) - t387 * r_i_i_C(2);
	t402 = -t404 - t427;
	t366 = t374 * t384 + t388 * t421;
	t400 = t374 * t388 - t384 * t421;
	t399 = t382 * t384 + t388 * t420;
	t398 = -t382 * t388 + t384 * t420;
	t397 = t412 - t415;
	t395 = qJD(5) * t404;
	t391 = qJD(3) + t384 * t395 + (-t425 * t384 + t403 * t388) * qJD(4);
	t363 = t399 * qJD(4) - t384 * t407;
	t362 = -qJD(1) * t412 - t386 * t413 + (qJD(2) * t382 + qJD(1)) * t415;
	t360 = t373 * qJD(1) + t374 * qJD(2);
	t359 = t372 * qJD(1) + t397 * qJD(2);
	t357 = t400 * qJD(4) - t359 * t384 + t388 * t409;
	t356 = t366 * qJD(4) + t359 * t388 + t384 * t409;
	t351 = t357 * t387 - t360 * t383 + (-t366 * t383 - t387 * t397) * qJD(5);
	t350 = -t357 * t383 - t360 * t387 + (-t366 * t387 + t383 * t397) * qJD(5);
	t1 = [-t361 * qJ(3) - t372 * qJD(3) - t403 * t354 - t402 * t362 + t425 * t428 + (t430 * r_i_i_C(1) - t429 * r_i_i_C(2)) * qJD(5) + (-t390 * pkin(1) - t426 * t421) * qJD(1), t402 * t359 - t392 * t360 - t374 * t396 - t391 * t397, -t359, -t403 * t356 - t425 * t357 + t400 * t395, t350 * r_i_i_C(1) - t351 * r_i_i_C(2); t357 * pkin(4) + t351 * r_i_i_C(1) + t350 * r_i_i_C(2) - t359 * qJ(3) + t374 * qJD(3) + t427 * t360 - t425 * t356 + (-pkin(1) * t386 + t426 * t419) * qJD(1), -t402 * t361 + t392 * t362 - t372 * t396 + t391 * t373, t361, -t425 * t354 + t401 * t395 - t403 * t428, (-t354 * t383 + t362 * t387) * r_i_i_C(1) + (-t354 * t387 - t362 * t383) * r_i_i_C(2) + (t429 * r_i_i_C(1) + t430 * r_i_i_C(2)) * qJD(5); 0, ((t392 * qJD(2) + t396) * t389 + (-t402 * qJD(2) + t391) * t385) * t381, t407, t425 * t363 - t399 * t395 + t403 * (t398 * qJD(4) + t388 * t407), (t363 * t383 + t387 * t408) * r_i_i_C(1) + (t363 * t387 - t383 * t408) * r_i_i_C(2) + ((-t383 * t422 + t387 * t398) * r_i_i_C(1) + (-t383 * t398 - t387 * t422) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end