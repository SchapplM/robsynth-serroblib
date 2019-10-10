% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(11));
	t50 = sin(pkin(11));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:30
	% EndTime: 2019-10-09 22:31:30
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(11));
	t182 = cos(pkin(11));
	t185 = sin(qJ(2));
	t183 = cos(pkin(6));
	t187 = cos(qJ(2));
	t195 = t183 * t187;
	t201 = -t180 * t185 + t182 * t195;
	t200 = pkin(8) + r_i_i_C(3);
	t181 = sin(pkin(6));
	t184 = sin(qJ(3));
	t198 = t181 * t184;
	t186 = cos(qJ(3));
	t197 = t181 * t186;
	t196 = t183 * t185;
	t193 = t184 * r_i_i_C(1) + t186 * r_i_i_C(2);
	t192 = t186 * r_i_i_C(1) - t184 * r_i_i_C(2) + pkin(2);
	t176 = t180 * t187 + t182 * t196;
	t191 = t180 * t195 + t182 * t185;
	t190 = t180 * t196 - t182 * t187;
	t189 = qJD(3) * t193;
	t188 = qJD(2) * t192;
	t173 = t191 * qJD(2);
	t171 = t201 * qJD(2);
	t1 = [0, -t200 * t173 + t190 * t188 + t191 * t189, t193 * t173 + ((-t180 * t198 + t186 * t190) * r_i_i_C(1) + (-t180 * t197 - t184 * t190) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, t200 * t171 - t176 * t188 - t201 * t189, -t193 * t171 + ((-t176 * t186 + t182 * t198) * r_i_i_C(1) + (t176 * t184 + t182 * t197) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (-t187 * t189 + (-t192 * t185 + t200 * t187) * qJD(2)) * t181, -t193 * t187 * t181 * qJD(2) + ((-t183 * t184 - t185 * t197) * r_i_i_C(1) + (-t183 * t186 + t185 * t198) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:30
	% EndTime: 2019-10-09 22:31:30
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (124->35), mult. (404->66), div. (0->0), fcn. (384->8), ass. (0->32)
	t231 = sin(pkin(11));
	t233 = cos(pkin(11));
	t236 = sin(qJ(2));
	t234 = cos(pkin(6));
	t238 = cos(qJ(2));
	t248 = t234 * t238;
	t259 = -t231 * t236 + t233 * t248;
	t249 = t234 * t236;
	t226 = t231 * t238 + t233 * t249;
	t237 = cos(qJ(3));
	t232 = sin(pkin(6));
	t235 = sin(qJ(3));
	t251 = t232 * t235;
	t258 = -t226 * t237 + t233 * t251;
	t254 = r_i_i_C(3) + qJ(4);
	t256 = pkin(3) + r_i_i_C(1);
	t257 = t254 * t235 + t256 * t237 + pkin(2);
	t255 = pkin(8) + r_i_i_C(2);
	t250 = t232 * t237;
	t245 = qJD(2) * t232 * t238;
	t242 = t231 * t249 - t233 * t238;
	t244 = t231 * t251 - t237 * t242;
	t243 = t231 * t248 + t233 * t236;
	t241 = t234 * t235 + t236 * t250;
	t240 = qJD(2) * t257;
	t239 = qJD(4) * t235 + (-t256 * t235 + t254 * t237) * qJD(3);
	t223 = t243 * qJD(2);
	t221 = t259 * qJD(2);
	t219 = t241 * qJD(3) + t235 * t245;
	t217 = t244 * qJD(3) - t223 * t235;
	t215 = -t258 * qJD(3) + t221 * t235;
	t1 = [0, -t255 * t223 - t239 * t243 + t242 * t240, t244 * qJD(4) + t254 * (-t223 * t237 + (t231 * t250 + t235 * t242) * qJD(3)) - t256 * t217, t217, 0, 0; 0, t255 * t221 - t226 * t240 + t239 * t259, -t258 * qJD(4) + t254 * (t221 * t237 + (-t226 * t235 - t233 * t250) * qJD(3)) - t256 * t215, t215, 0, 0; 0, (t239 * t238 + (-t257 * t236 + t255 * t238) * qJD(2)) * t232, t241 * qJD(4) + t254 * (t237 * t245 + (t234 * t237 - t236 * t251) * qJD(3)) - t256 * t219, t219, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:30
	% EndTime: 2019-10-09 22:31:31
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (331->67), mult. (1044->117), div. (0->0), fcn. (1080->10), ass. (0->44)
	t260 = sin(pkin(11));
	t262 = cos(pkin(11));
	t269 = cos(qJ(2));
	t263 = cos(pkin(6));
	t266 = sin(qJ(2));
	t289 = t263 * t266;
	t251 = t260 * t269 + t262 * t289;
	t265 = sin(qJ(3));
	t261 = sin(pkin(6));
	t268 = cos(qJ(3));
	t290 = t261 * t268;
	t240 = t251 * t265 + t262 * t290;
	t291 = t261 * t265;
	t241 = t251 * t268 - t262 * t291;
	t264 = sin(qJ(5));
	t267 = cos(qJ(5));
	t299 = ((t240 * t264 + t241 * t267) * r_i_i_C(1) + (t240 * t267 - t241 * t264) * r_i_i_C(2)) * qJD(5);
	t272 = t260 * t289 - t262 * t269;
	t243 = t260 * t291 - t268 * t272;
	t274 = t260 * t290 + t265 * t272;
	t298 = ((t243 * t267 - t264 * t274) * r_i_i_C(1) + (-t243 * t264 - t267 * t274) * r_i_i_C(2)) * qJD(5);
	t254 = -t263 * t268 + t266 * t291;
	t255 = t263 * t265 + t266 * t290;
	t297 = ((t254 * t264 + t255 * t267) * r_i_i_C(1) + (t254 * t267 - t255 * t264) * r_i_i_C(2)) * qJD(5);
	t275 = t267 * r_i_i_C(1) - t264 * r_i_i_C(2) + pkin(3) + pkin(4);
	t276 = t264 * r_i_i_C(1) + t267 * r_i_i_C(2) + qJ(4);
	t293 = t276 * t265 + t275 * t268 + pkin(2);
	t288 = t263 * t269;
	t287 = qJD(2) * t266;
	t286 = r_i_i_C(3) + pkin(9) - pkin(8);
	t284 = t262 * t288;
	t283 = qJD(2) * t261 * t269;
	t273 = t260 * t288 + t262 * t266;
	t271 = qJD(2) * t293;
	t270 = t265 * qJD(4) + ((-t264 * t268 + t265 * t267) * r_i_i_C(1) + (-t264 * t265 - t267 * t268) * r_i_i_C(2)) * qJD(5) + (-t275 * t265 + t276 * t268) * qJD(3);
	t248 = t273 * qJD(2);
	t246 = -qJD(2) * t284 + t260 * t287;
	t245 = -t254 * qJD(3) + t268 * t283;
	t244 = t255 * qJD(3) + t265 * t283;
	t239 = t274 * qJD(3) - t248 * t268;
	t238 = t243 * qJD(3) - t248 * t265;
	t237 = -t240 * qJD(3) - t246 * t268;
	t236 = t241 * qJD(3) - t246 * t265;
	t1 = [0, t286 * t248 - t270 * t273 + t272 * t271, t243 * qJD(4) - t275 * t238 + t276 * t239 + t298, t238, (t238 * t267 - t239 * t264) * r_i_i_C(1) + (-t238 * t264 - t239 * t267) * r_i_i_C(2) - t298, 0; 0, t286 * t246 - t251 * t271 + t270 * (-t260 * t266 + t284), t241 * qJD(4) - t275 * t236 + t276 * t237 + t299, t236, (t236 * t267 - t237 * t264) * r_i_i_C(1) + (-t236 * t264 - t237 * t267) * r_i_i_C(2) - t299, 0; 0, (-t293 * t287 + (-t286 * qJD(2) + t270) * t269) * t261, t255 * qJD(4) - t275 * t244 + t276 * t245 + t297, t244, (t244 * t267 - t245 * t264) * r_i_i_C(1) + (-t244 * t264 - t245 * t267) * r_i_i_C(2) - t297, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:34
	% EndTime: 2019-10-09 22:31:35
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (905->115), mult. (2786->204), div. (0->0), fcn. (3016->12), ass. (0->72)
	t517 = sin(qJ(5));
	t518 = sin(qJ(3));
	t521 = cos(qJ(5));
	t522 = cos(qJ(3));
	t534 = t517 * t518 + t521 * t522;
	t565 = qJD(5) - qJD(3);
	t524 = t565 * t534;
	t515 = cos(pkin(6));
	t519 = sin(qJ(2));
	t513 = sin(pkin(6));
	t557 = t513 * t522;
	t506 = t515 * t518 + t519 * t557;
	t523 = cos(qJ(2));
	t556 = t513 * t523;
	t548 = qJD(2) * t556;
	t494 = t506 * qJD(3) + t518 * t548;
	t558 = t513 * t518;
	t505 = -t515 * t522 + t519 * t558;
	t495 = -t505 * qJD(3) + t522 * t548;
	t536 = t505 * t521 - t506 * t517;
	t472 = t536 * qJD(5) + t494 * t517 + t495 * t521;
	t489 = t505 * t517 + t506 * t521;
	t516 = sin(qJ(6));
	t520 = cos(qJ(6));
	t545 = -t516 * r_i_i_C(1) - t520 * r_i_i_C(2);
	t528 = qJD(6) * t545;
	t546 = -t520 * r_i_i_C(1) + t516 * r_i_i_C(2);
	t533 = pkin(5) - t546;
	t560 = r_i_i_C(3) + pkin(10);
	t574 = -t533 * (t489 * qJD(5) - t494 * t521 + t495 * t517) + t560 * t472 + t536 * t528;
	t512 = sin(pkin(11));
	t514 = cos(pkin(11));
	t555 = t515 * t519;
	t529 = t512 * t555 - t514 * t523;
	t493 = t512 * t558 - t522 * t529;
	t554 = t515 * t523;
	t530 = t512 * t554 + t514 * t519;
	t499 = t530 * qJD(2);
	t483 = t493 * qJD(3) - t499 * t518;
	t531 = t512 * t557 + t518 * t529;
	t484 = t531 * qJD(3) - t499 * t522;
	t539 = -t493 * t517 - t521 * t531;
	t464 = t539 * qJD(5) + t483 * t517 + t484 * t521;
	t480 = t493 * t521 - t517 * t531;
	t573 = -t533 * (t480 * qJD(5) - t483 * t521 + t484 * t517) + t560 * t464 + t539 * t528;
	t502 = t512 * t523 + t514 * t555;
	t491 = t502 * t522 - t514 * t558;
	t550 = t514 * t554;
	t553 = qJD(2) * t519;
	t497 = -qJD(2) * t550 + t512 * t553;
	t481 = t491 * qJD(3) - t497 * t518;
	t490 = t502 * t518 + t514 * t557;
	t482 = -t490 * qJD(3) - t497 * t522;
	t540 = t490 * t521 - t491 * t517;
	t460 = t540 * qJD(5) + t481 * t517 + t482 * t521;
	t477 = t490 * t517 + t491 * t521;
	t572 = -t533 * (t477 * qJD(5) - t481 * t521 + t482 * t517) + t560 * t460 + t540 * t528;
	t561 = pkin(3) + pkin(4);
	t527 = t518 * qJ(4) + t561 * t522 + pkin(2);
	t552 = qJD(6) * t534 * t556;
	t549 = t513 * t553;
	t535 = t517 * t522 - t518 * t521;
	t532 = -pkin(8) + pkin(9) - t545;
	t526 = t518 * qJD(4) + (qJ(4) * t522 - t561 * t518) * qJD(3);
	t525 = t565 * t535;
	t501 = -t512 * t519 + t550;
	t500 = t529 * qJD(2);
	t498 = t502 * qJD(2);
	t486 = t534 * t530;
	t485 = t534 * t501;
	t474 = (-t525 * t523 - t534 * t553) * t513;
	t1 = [0, t533 * (t534 * t500 + t525 * t530) + t532 * t499 + t560 * (t535 * t500 - t524 * t530) + ((t486 * t516 + t520 * t529) * r_i_i_C(1) + (t486 * t520 - t516 * t529) * r_i_i_C(2)) * qJD(6) - t526 * t530 + t527 * t500, t484 * qJ(4) + t493 * qJD(4) - t561 * t483 - t573, t483, t573, (-t464 * t516 + t500 * t520) * r_i_i_C(1) + (-t464 * t520 - t500 * t516) * r_i_i_C(2) + ((-t480 * t520 + t516 * t530) * r_i_i_C(1) + (t480 * t516 + t520 * t530) * r_i_i_C(2)) * qJD(6); 0, t533 * (-t534 * t498 - t525 * t501) + t532 * t497 + t560 * (-t535 * t498 + t524 * t501) + ((-t485 * t516 - t502 * t520) * r_i_i_C(1) + (-t485 * t520 + t502 * t516) * r_i_i_C(2)) * qJD(6) + t526 * t501 - t527 * t498, t482 * qJ(4) + t491 * qJD(4) - t561 * t481 - t572, t481, t572, (-t460 * t516 - t498 * t520) * r_i_i_C(1) + (-t460 * t520 + t498 * t516) * r_i_i_C(2) + ((-t477 * t520 - t501 * t516) * r_i_i_C(1) + (t477 * t516 - t501 * t520) * r_i_i_C(2)) * qJD(6); 0, (t474 * t520 - t516 * t552) * r_i_i_C(1) + (-t474 * t516 - t520 * t552) * r_i_i_C(2) + t474 * pkin(5) + (-t560 * (-t523 * t524 + t535 * t553) + t546 * t519 * qJD(6) + t526 * t523 + (-t527 * t519 - t532 * t523) * qJD(2)) * t513, t495 * qJ(4) + t506 * qJD(4) - t561 * t494 - t574, t494, t574, (-t472 * t516 - t520 * t549) * r_i_i_C(1) + (-t472 * t520 + t516 * t549) * r_i_i_C(2) + ((-t489 * t520 - t516 * t556) * r_i_i_C(1) + (t489 * t516 - t520 * t556) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end