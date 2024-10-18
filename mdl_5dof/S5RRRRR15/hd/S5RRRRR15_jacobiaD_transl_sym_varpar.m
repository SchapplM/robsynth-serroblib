% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR15
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR15_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR15_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR15_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(5));
	t151 = t136 * (pkin(8) + r_i_i_C(3));
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
	t135 = qJD(1) * t142 + qJD(2) * t145;
	t134 = qJD(1) * t144 + qJD(2) * t143;
	t133 = qJD(1) * t143 + qJD(2) * t144;
	t132 = qJD(1) * t145 + qJD(2) * t142;
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:08
	% EndTime: 2024-09-27 22:28:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (277->47), mult. (231->76), div. (0->0), fcn. (176->12), ass. (0->49)
	t227 = sin(pkin(5));
	t229 = sin(qJ(2));
	t253 = pkin(2) * qJD(2);
	t237 = t229 * t253;
	t258 = qJD(1) * t227 * r_i_i_C(3) - t237;
	t226 = qJ(2) + qJ(3);
	t222 = pkin(5) - t226;
	t219 = cos(t222);
	t225 = qJD(2) + qJD(3);
	t215 = t225 * t219;
	t221 = pkin(5) + t226;
	t218 = cos(t221);
	t252 = t225 * t218;
	t209 = t215 + t252;
	t216 = sin(t221);
	t217 = sin(t222);
	t212 = t216 - t217;
	t232 = cos(qJ(1));
	t223 = sin(t226);
	t245 = t232 * t223;
	t224 = cos(t226);
	t230 = sin(qJ(1));
	t248 = t230 * t224;
	t256 = -t230 / 0.2e1;
	t257 = (t248 + t232 * t212 / 0.2e1) * qJD(1) - t209 * t256 + t225 * t245;
	t255 = t230 / 0.2e1;
	t254 = -t232 / 0.2e1;
	t251 = t229 * t230;
	t250 = t229 * t232;
	t249 = t230 * t223;
	t231 = cos(qJ(2));
	t247 = t230 * t231;
	t246 = t231 * t232;
	t244 = t232 * t224;
	t210 = t212 * t225;
	t213 = t219 + t218;
	t205 = -t225 * t244 - t210 * t256 + (t213 * t254 + t249) * qJD(1);
	t243 = t205 * r_i_i_C(1) + t257 * r_i_i_C(2);
	t206 = t225 * t248 - t210 * t254 + (t213 * t255 + t245) * qJD(1);
	t233 = t225 * t249 + t209 * t254 + (t212 * t255 - t244) * qJD(1);
	t242 = -t206 * r_i_i_C(1) + t233 * r_i_i_C(2);
	t241 = (t252 / 0.2e1 - t215 / 0.2e1) * r_i_i_C(1) + (-t216 / 0.2e1 - t217 / 0.2e1) * t225 * r_i_i_C(2);
	t240 = qJD(1) * t230;
	t239 = qJD(1) * t232;
	t228 = cos(pkin(5));
	t236 = t228 * t231 * t253;
	t220 = t231 * pkin(2) + pkin(1);
	t211 = t228 * t229 * pkin(2) - t227 * (pkin(8) + pkin(9));
	t1 = [t233 * r_i_i_C(1) + t206 * r_i_i_C(2) + t211 * t240 - t220 * t239 - t258 * t230 - t232 * t236, ((t228 * t251 - t246) * qJD(2) + (-t228 * t246 + t251) * qJD(1)) * pkin(2) + t243, t243, 0, 0; -t257 * r_i_i_C(1) + t205 * r_i_i_C(2) - t211 * t239 - t220 * t240 - t230 * t236 + t258 * t232, ((-t228 * t250 - t247) * qJD(2) + (-t228 * t247 - t250) * qJD(1)) * pkin(2) + t242, t242, 0, 0; 0, -t227 * t237 + t241, t241, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:08
	% EndTime: 2024-09-27 22:28:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (610->75), mult. (492->117), div. (0->0), fcn. (354->16), ass. (0->63)
	t240 = qJ(2) + qJ(3);
	t238 = qJ(4) + t240;
	t230 = pkin(5) + t238;
	t275 = cos(t230) / 0.2e1;
	t239 = qJD(2) + qJD(3);
	t235 = qJD(4) + t239;
	t274 = -t235 / 0.2e1;
	t231 = pkin(5) - t238;
	t273 = cos(t231);
	t272 = sin(t230);
	t236 = sin(t240);
	t271 = pkin(3) * t236;
	t228 = sin(t231);
	t270 = t235 * t228;
	t245 = sin(qJ(1));
	t269 = t235 * t245;
	t248 = cos(qJ(1));
	t268 = t235 * t248;
	t237 = cos(t240);
	t267 = t237 * t239;
	t243 = sin(qJ(3));
	t244 = sin(qJ(2));
	t266 = t243 * t244;
	t247 = cos(qJ(2));
	t265 = t243 * t247;
	t220 = t272 * t274;
	t213 = t270 / 0.2e1 + t220;
	t217 = t273 / 0.2e1 + t275;
	t232 = sin(t238);
	t233 = cos(t238);
	t204 = -t233 * t268 - t245 * t213 + (-t217 * t248 + t232 * t245) * qJD(1);
	t216 = t272 / 0.2e1 - t228 / 0.2e1;
	t212 = t217 * t235;
	t252 = -t245 * t212 - t232 * t268;
	t264 = t204 * r_i_i_C(1) + ((t216 * t248 + t233 * t245) * qJD(1) - t252) * r_i_i_C(2);
	t205 = t233 * t269 - t248 * t213 + (t217 * t245 + t232 * t248) * qJD(1);
	t256 = -t248 * t212 + t232 * t269;
	t263 = -t205 * r_i_i_C(1) + ((t216 * t245 - t233 * t248) * qJD(1) + t256) * r_i_i_C(2);
	t262 = (t235 * t275 + t273 * t274) * r_i_i_C(1) + (t220 - t270 / 0.2e1) * r_i_i_C(2);
	t261 = qJD(1) * t245;
	t260 = qJD(1) * t248;
	t259 = qJD(2) * t244;
	t258 = qJD(2) * t247;
	t257 = -t247 * pkin(2) - pkin(3) * t237 - t233 * r_i_i_C(1) - pkin(1);
	t246 = cos(qJ(3));
	t255 = -t244 * t246 - t265;
	t254 = t246 * t247 - t266;
	t234 = t246 * pkin(3) + pkin(2);
	t241 = sin(pkin(5));
	t242 = cos(pkin(5));
	t253 = t216 * r_i_i_C(1) + (pkin(3) * t265 + t244 * t234) * t242 + (-r_i_i_C(3) - pkin(8) - pkin(9) - pkin(10)) * t241;
	t251 = t255 * qJD(3);
	t250 = pkin(3) * (t255 * qJD(2) + t251);
	t249 = -t234 * t259 + (-t243 * t258 + t251) * pkin(3);
	t225 = -t244 * pkin(2) - t271;
	t219 = -pkin(2) * t258 - pkin(3) * t267;
	t218 = -pkin(2) * t259 - t239 * t271;
	t215 = t254 * t242 * pkin(3);
	t214 = (-pkin(3) * t266 + t234 * t247) * t242;
	t208 = t242 * t250;
	t207 = (t234 * t258 + (t254 * qJD(3) - t243 * t259) * pkin(3)) * t242;
	t206 = t249 * t242;
	t1 = [t256 * r_i_i_C(1) + t205 * r_i_i_C(2) - t245 * t218 - t248 * t207 + (t253 * t245 + t257 * t248) * qJD(1), -t245 * t206 + t248 * t219 + (-t214 * t248 - t225 * t245) * qJD(1) + t264, -t215 * t260 - t245 * t208 + (t236 * t261 - t248 * t267) * pkin(3) + t264, t264, 0; t252 * r_i_i_C(1) + t204 * r_i_i_C(2) + t248 * t218 - t245 * t207 + (t257 * t245 - t253 * t248) * qJD(1), t248 * t206 + t245 * t219 + (-t214 * t245 + t225 * t248) * qJD(1) + t263, -t215 * t261 + t248 * t208 + (-t236 * t260 - t245 * t267) * pkin(3) + t263, t263, 0; 0, t249 * t241 + t262, t241 * t250 + t262, t262, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:09
	% EndTime: 2024-09-27 22:28:09
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (1448->159), mult. (2762->261), div. (0->0), fcn. (2481->16), ass. (0->121)
	t483 = cos(qJ(4));
	t474 = sin(pkin(6));
	t552 = pkin(11) * t474;
	t553 = sin(qJ(4));
	t466 = -pkin(4) * t483 - t553 * t552;
	t465 = pkin(3) - t466;
	t479 = sin(qJ(3));
	t527 = qJD(3) * t479;
	t467 = -t553 * pkin(4) + t483 * t552;
	t484 = cos(qJ(3));
	t464 = t467 * t484;
	t488 = t466 * qJD(4);
	t497 = t467 * qJD(4);
	t555 = -qJD(3) * t464 - t479 * t488 - t497 * t484;
	t427 = t465 * t527 + t555;
	t542 = t467 * t479;
	t447 = -t465 * t484 - t542;
	t556 = -t479 * t497 + t484 * t488;
	t428 = t447 * qJD(3) + t556;
	t480 = sin(qJ(2));
	t485 = cos(qJ(2));
	t559 = t427 * t480 + t428 * t485;
	t446 = pkin(2) - t447;
	t448 = -t479 * t465 + t464;
	t543 = t448 * t480;
	t501 = t446 * t485 + t543;
	t408 = -t501 * qJD(2) + t559;
	t445 = t448 * t485;
	t426 = t447 * t480 + t445;
	t558 = -t446 * t480 + t445;
	t473 = qJ(2) + qJ(3) + qJ(4);
	t470 = cos(t473);
	t475 = sin(pkin(5));
	t476 = cos(pkin(6));
	t477 = cos(pkin(5));
	t540 = t474 * t477;
	t557 = -(t470 * t540 + t475 * t476) * r_i_i_C(3) - t475 * (t476 * pkin(11) + pkin(8) + pkin(9) + pkin(10)) - t558 * t477;
	t449 = t466 * t484 - t542;
	t429 = t449 * qJD(3) + t556;
	t430 = -t466 * t527 + t555;
	t450 = t466 * t479 + t464;
	t431 = t449 * t480 + t450 * t485;
	t554 = t431 * qJD(2) + t429 * t480 - t430 * t485;
	t478 = sin(qJ(5));
	t551 = r_i_i_C(1) * t478;
	t482 = cos(qJ(5));
	t550 = r_i_i_C(2) * t482;
	t549 = r_i_i_C(3) * t474;
	t547 = t427 * t485;
	t546 = t428 * t480;
	t541 = t474 * t475;
	t539 = t476 * t478;
	t538 = t476 * t482;
	t481 = sin(qJ(1));
	t537 = t477 * t481;
	t486 = cos(qJ(1));
	t536 = t477 * t486;
	t535 = t481 * t478;
	t534 = t482 * t481;
	t533 = t486 * t478;
	t532 = t486 * t482;
	t530 = qJD(1) * t481;
	t529 = qJD(1) * t486;
	t528 = qJD(2) * t480;
	t526 = qJD(5) * t481;
	t525 = qJD(5) * t486;
	t522 = t477 * t533;
	t521 = t477 * t535;
	t520 = t477 * t534;
	t519 = t477 * t532;
	t469 = sin(t473);
	t472 = qJD(2) + qJD(3) + qJD(4);
	t492 = -t469 * t536 - t470 * t481;
	t456 = t476 * t534 + t522;
	t457 = t476 * t533 + t520;
	t489 = t476 * t520 + t533;
	t507 = t456 * qJD(1) + t457 * qJD(5) + t472 * t489;
	t455 = t476 * t535 - t519;
	t458 = -t476 * t532 + t521;
	t441 = t455 * qJD(1) + t458 * qJD(5);
	t451 = t476 * t521 - t532;
	t508 = t451 * t472 + t441;
	t453 = t476 * t522 + t534;
	t439 = t453 * qJD(1) + t489 * qJD(5);
	t510 = t457 * t472 + t439;
	t454 = t476 * t519 - t535;
	t513 = -t454 * qJD(1) + t451 * qJD(5) + t458 * t472;
	t518 = (-t513 * t469 + t507 * t470) * r_i_i_C(2) + (t510 * t469 + t508 * t470) * r_i_i_C(1) + ((-t469 * t486 - t470 * t537) * t472 + t492 * qJD(1)) * t549;
	t491 = -t469 * t537 + t470 * t486;
	t506 = t457 * qJD(1) + t456 * qJD(5) + t453 * t472;
	t440 = t458 * qJD(1) + t455 * qJD(5);
	t509 = -t454 * t472 + t440;
	t438 = t489 * qJD(1) + t453 * qJD(5);
	t511 = t456 * t472 + t438;
	t512 = t451 * qJD(1) - t454 * qJD(5) + t455 * t472;
	t517 = (t512 * t469 - t506 * t470) * r_i_i_C(1) + (t511 * t469 + t509 * t470) * r_i_i_C(2) + ((-t469 * t481 + t470 * t536) * t472 + t491 * qJD(1)) * t549;
	t493 = t469 * t539 - t470 * t482;
	t494 = -t469 * t538 - t470 * t478;
	t495 = -t469 * t482 - t470 * t539;
	t496 = t469 * t478 - t470 * t538;
	t516 = t472 * t470 * r_i_i_C(3) * t541 + ((t493 * qJD(5) + t496 * t472) * r_i_i_C(2) + (t494 * qJD(5) + t495 * t472) * r_i_i_C(1)) * t475;
	t515 = t482 * t526;
	t514 = t478 * t525;
	t505 = -t469 * t549 - pkin(1) - t501;
	t504 = t516 + (qJD(2) * t445 + t546) * t475;
	t503 = t546 - t547;
	t499 = t447 * t485 - t543;
	t498 = t449 * t485 - t450 * t480;
	t409 = qJD(2) * t558 + t503;
	t422 = t498 * t477;
	t420 = t499 * t477;
	t419 = t501 * t477;
	t413 = (-t478 * t530 + t482 * t525) * t541 + t512 * t470 + t506 * t469;
	t412 = (-t478 * t526 + t482 * t529) * t541 + t513 * t470 + t507 * t469;
	t411 = t498 * qJD(2) + t429 * t485 + t430 * t480;
	t410 = t499 * qJD(2) + t559;
	t407 = t554 * t477;
	t406 = (-t426 * qJD(2) - t503) * t477;
	t405 = t408 * t477;
	t404 = t409 * t477;
	t1 = [t413 * r_i_i_C(1) + (t438 * t470 - t440 * t469 - t514 * t541) * r_i_i_C(2) + t405 * t486 - t409 * t481 + ((t454 * t469 + t456 * t470) * r_i_i_C(2) + t492 * t549) * t472 + (t505 * t486 + (-t541 * t550 + t557) * t481) * qJD(1), -t404 * t481 + t408 * t486 + (-t419 * t486 - t481 * t558) * qJD(1) + t518, t406 * t481 + t410 * t486 + (t420 * t486 - t426 * t481) * qJD(1) + t518, -t407 * t481 + t411 * t486 + (t422 * t486 - t431 * t481) * qJD(1) + t518, t412 * r_i_i_C(1) + ((-t478 * t529 - t515) * t541 + t510 * t470 - t508 * t469) * r_i_i_C(2); (-t439 * t470 + t441 * t469 + t515 * t541) * r_i_i_C(1) + t412 * r_i_i_C(2) + t405 * t481 + t409 * t486 + ((t451 * t469 - t457 * t470) * r_i_i_C(1) + t491 * t549) * t472 + (t505 * t481 + (t541 * t551 - t557) * t486) * qJD(1), t404 * t486 + t408 * t481 + (-t419 * t481 + t486 * t558) * qJD(1) + t517, -t406 * t486 + t410 * t481 + (t420 * t481 + t426 * t486) * qJD(1) + t517, t407 * t486 + t411 * t481 + (t422 * t481 + t431 * t486) * qJD(1) + t517, ((t482 * t530 + t514) * t541 - t511 * t470 + t509 * t469) * r_i_i_C(1) + t413 * r_i_i_C(2); 0, (-t446 * t528 - t547) * t475 + t504, (t447 * t528 - t547) * t475 + t504, t554 * t475 + t516, (-t550 - t551) * qJD(5) * t540 + ((t494 * r_i_i_C(1) + t493 * r_i_i_C(2)) * t472 + (t495 * r_i_i_C(1) + t496 * r_i_i_C(2)) * qJD(5)) * t475;];
	JaD_transl = t1;
end