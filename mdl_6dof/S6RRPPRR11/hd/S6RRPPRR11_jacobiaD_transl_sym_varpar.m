% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
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
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
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
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:52
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (137->38), mult. (400->56), div. (0->0), fcn. (372->8), ass. (0->29)
	t203 = sin(pkin(11));
	t204 = sin(pkin(6));
	t205 = cos(pkin(11));
	t223 = t204 * (r_i_i_C(1) * t205 - r_i_i_C(2) * t203 + pkin(3) + pkin(8));
	t207 = sin(qJ(2));
	t208 = sin(qJ(1));
	t222 = t208 * t207;
	t209 = cos(qJ(2));
	t221 = t208 * t209;
	t210 = cos(qJ(1));
	t220 = t210 * t207;
	t219 = t210 * t209;
	t218 = qJD(2) * t207;
	t217 = qJD(2) * t209;
	t216 = r_i_i_C(3) + qJ(4) + pkin(2);
	t206 = cos(pkin(6));
	t215 = t206 * t222;
	t214 = t206 * t219;
	t213 = qJD(2) * t206 + qJD(1);
	t212 = t203 * r_i_i_C(1) + t205 * r_i_i_C(2) + qJ(3);
	t197 = t206 * t221 + t220;
	t196 = t206 * t220 + t221;
	t198 = t215 - t219;
	t195 = -t214 + t222;
	t194 = -qJD(1) * t215 - t208 * t218 + t213 * t219;
	t193 = t197 * qJD(1) + t196 * qJD(2);
	t192 = t196 * qJD(1) + t197 * qJD(2);
	t191 = -qJD(1) * t214 - t210 * t217 + t213 * t222;
	t1 = [-t195 * qJD(3) - t196 * qJD(4) - t216 * t194 - t212 * t193 + (-t210 * pkin(1) - t208 * t223) * qJD(1), -t198 * qJD(3) - t197 * qJD(4) + t216 * t191 - t212 * t192, -t191, -t192, 0, 0; t197 * qJD(3) - t198 * qJD(4) - t216 * t192 - t212 * t191 + (-t208 * pkin(1) + t210 * t223) * qJD(1), t196 * qJD(3) - t195 * qJD(4) - t216 * t193 + t212 * t194, t193, t194, 0, 0; 0, (qJD(3) * t207 + qJD(4) * t209 + (-t216 * t207 + t212 * t209) * qJD(2)) * t204, t204 * t218, t204 * t217, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:52
	% EndTime: 2019-10-10 09:53:52
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (268->63), mult. (628->102), div. (0->0), fcn. (595->10), ass. (0->44)
	t252 = pkin(11) + qJ(5);
	t250 = sin(t252);
	t251 = cos(t252);
	t264 = r_i_i_C(1) * t251 - r_i_i_C(2) * t250;
	t261 = t264 * qJD(5) + qJD(3);
	t284 = pkin(8) + cos(pkin(11)) * pkin(4) + pkin(3);
	t254 = sin(pkin(6));
	t258 = sin(qJ(1));
	t283 = t254 * t258;
	t259 = cos(qJ(2));
	t282 = t254 * t259;
	t260 = cos(qJ(1));
	t281 = t254 * t260;
	t257 = sin(qJ(2));
	t280 = t258 * t257;
	t279 = t258 * t259;
	t278 = t260 * t257;
	t277 = t260 * t259;
	t276 = qJD(1) * t258;
	t275 = qJD(1) * t260;
	t274 = qJD(2) * t257;
	t273 = qJD(2) * t259;
	t272 = -r_i_i_C(3) - pkin(9) - qJ(4) - pkin(2);
	t255 = cos(pkin(6));
	t271 = t255 * t280;
	t270 = t255 * t277;
	t269 = t254 * t276;
	t268 = t254 * t275;
	t267 = t254 * t274;
	t266 = -sin(pkin(11)) * pkin(4) - qJ(3);
	t265 = qJD(2) * t255 + qJD(1);
	t263 = -t250 * r_i_i_C(1) - t251 * r_i_i_C(2);
	t241 = t255 * t279 + t278;
	t240 = t255 * t278 + t279;
	t262 = -t263 - t266;
	t242 = -t271 + t277;
	t239 = -t270 + t280;
	t238 = -qJD(1) * t271 - t258 * t274 + t265 * t277;
	t237 = t241 * qJD(1) + t240 * qJD(2);
	t236 = t240 * qJD(1) + t241 * qJD(2);
	t235 = -qJD(1) * t270 - t260 * t273 + t265 * t280;
	t234 = t251 * t268 - t235 * t250 + (t241 * t251 - t250 * t283) * qJD(5);
	t233 = -t250 * t268 - t235 * t251 + (-t241 * t250 - t251 * t283) * qJD(5);
	t1 = [-pkin(1) * t275 - t240 * qJD(4) - t261 * t239 - t262 * t237 + t272 * t238 + (t263 * t260 * qJD(5) + (-t264 - t284) * t276) * t254, -t241 * qJD(4) - t272 * t235 - t262 * t236 + t261 * t242, -t235, -t236, t233 * r_i_i_C(1) - t234 * r_i_i_C(2), 0; t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t241 * qJD(3) + t242 * qJD(4) + t266 * t235 + t272 * t236 + (-pkin(1) * t258 + t284 * t281) * qJD(1), -t239 * qJD(4) + t272 * t237 + t262 * t238 + t261 * t240, t237, t238, (t237 * t251 - t250 * t269) * r_i_i_C(1) + (-t237 * t250 - t251 * t269) * r_i_i_C(2) + ((-t239 * t250 + t251 * t281) * r_i_i_C(1) + (-t239 * t251 - t250 * t281) * r_i_i_C(2)) * qJD(5), 0; 0, (qJD(4) * t259 + t261 * t257 + (t272 * t257 + t262 * t259) * qJD(2)) * t254, t267, t254 * t273, t264 * t267 + ((t250 * t282 - t251 * t255) * r_i_i_C(1) + (t250 * t255 + t251 * t282) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:54
	% EndTime: 2019-10-10 09:53:55
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (702->100), mult. (1518->156), div. (0->0), fcn. (1517->12), ass. (0->66)
	t402 = cos(pkin(6));
	t408 = cos(qJ(2));
	t409 = cos(qJ(1));
	t437 = t409 * t408;
	t431 = t402 * t437;
	t405 = sin(qJ(2));
	t406 = sin(qJ(1));
	t440 = t406 * t405;
	t384 = -t431 + t440;
	t399 = pkin(11) + qJ(5);
	t397 = sin(t399);
	t398 = cos(t399);
	t401 = sin(pkin(6));
	t441 = t401 * t409;
	t432 = t398 * t441;
	t380 = -t384 * t397 + t432;
	t438 = t409 * t405;
	t439 = t406 * t408;
	t385 = t402 * t438 + t439;
	t404 = sin(qJ(6));
	t407 = cos(qJ(6));
	t452 = -t380 * t404 - t385 * t407;
	t451 = t380 * t407 - t385 * t404;
	t386 = t402 * t439 + t438;
	t373 = t386 * qJD(1) + t385 * qJD(2);
	t436 = qJD(1) * t401;
	t430 = t406 * t436;
	t450 = (qJD(5) * t384 + t430) * t397 - qJD(5) * t432 - t373 * t398;
	t420 = t384 * t398 + t397 * t441;
	t366 = t420 * qJD(5) + t373 * t397 + t398 * t430;
	t423 = t407 * r_i_i_C(1) - t404 * r_i_i_C(2);
	t412 = -t423 * qJD(6) - qJD(4);
	t421 = pkin(5) + t423;
	t426 = -sin(pkin(11)) * pkin(4) - qJ(3);
	t449 = r_i_i_C(3) + pkin(10);
	t411 = t421 * t397 - t449 * t398 - t426;
	t448 = -pkin(2) - pkin(9) - qJ(4);
	t447 = pkin(8) + cos(pkin(11)) * pkin(4) + pkin(3);
	t444 = t401 * t405;
	t443 = t401 * t406;
	t442 = t401 * t408;
	t435 = qJD(2) * t405;
	t434 = qJD(2) * t408;
	t433 = t402 * t440;
	t429 = t409 * t436;
	t428 = t401 * t435;
	t427 = t401 * t434;
	t424 = qJD(2) * t402 + qJD(1);
	t422 = -t404 * r_i_i_C(1) - t407 * r_i_i_C(2);
	t378 = t386 * t397 + t398 * t443;
	t419 = t386 * t398 - t397 * t443;
	t418 = t402 * t397 + t398 * t442;
	t417 = t397 * t442 - t402 * t398;
	t416 = qJD(6) * t422;
	t415 = -t422 - t448;
	t410 = qJD(3) + t397 * t416 + (t449 * t397 + t421 * t398) * qJD(5);
	t387 = -t433 + t437;
	t375 = t418 * qJD(5) - t397 * t428;
	t374 = -qJD(1) * t433 - t406 * t435 + t424 * t437;
	t372 = t385 * qJD(1) + t386 * qJD(2);
	t371 = -qJD(1) * t431 - t409 * t434 + t424 * t440;
	t369 = t419 * qJD(5) - t371 * t397 + t398 * t429;
	t368 = t378 * qJD(5) + t371 * t398 + t397 * t429;
	t363 = t369 * t407 - t372 * t404 + (-t378 * t404 + t387 * t407) * qJD(6);
	t362 = -t369 * t404 - t372 * t407 + (-t378 * t407 - t387 * t404) * qJD(6);
	t1 = [-t384 * qJD(3) - t385 * qJD(4) - t421 * t366 - t415 * t374 + t426 * t373 - t449 * t450 + (t452 * r_i_i_C(1) - t451 * r_i_i_C(2)) * qJD(6) + (-t409 * pkin(1) - t447 * t443) * qJD(1), t415 * t371 - t411 * t372 + t412 * t386 + t410 * t387, -t371, -t372, -t421 * t368 + t449 * t369 + t419 * t416, t362 * r_i_i_C(1) - t363 * r_i_i_C(2); t369 * pkin(5) + t363 * r_i_i_C(1) + t362 * r_i_i_C(2) + t386 * qJD(3) + t387 * qJD(4) + t448 * t372 + t426 * t371 + t449 * t368 + (-pkin(1) * t406 + t447 * t441) * qJD(1), -t415 * t373 + t411 * t374 + t412 * t384 + t410 * t385, t373, t374, t449 * t366 + t420 * t416 - t421 * t450, (-t366 * t404 + t374 * t407) * r_i_i_C(1) + (-t366 * t407 - t374 * t404) * r_i_i_C(2) + (t451 * r_i_i_C(1) + t452 * r_i_i_C(2)) * qJD(6); 0, ((t411 * qJD(2) - t412) * t408 + (-t415 * qJD(2) + t410) * t405) * t401, t428, t427, -t449 * t375 - t418 * t416 + t421 * (t417 * qJD(5) + t398 * t428), (t375 * t404 + t407 * t427) * r_i_i_C(1) + (t375 * t407 - t404 * t427) * r_i_i_C(2) + ((-t404 * t444 + t407 * t417) * r_i_i_C(1) + (-t404 * t417 - t407 * t444) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end