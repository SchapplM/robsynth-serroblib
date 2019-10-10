% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:19
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(6));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(6));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0, 0; -t153, -t154, 0, 0, 0, 0; 0, -t158 * t166, 0, 0, 0, 0; t154, t153, 0, 0, 0, 0; t152, t155, 0, 0, 0, 0; 0, -t160 * t166, 0, 0, 0, 0; -t159 * t167, 0, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:19:00
	% EndTime: 2019-10-10 10:19:00
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (43->18), mult. (148->39), div. (0->0), fcn. (148->8), ass. (0->25)
	t231 = sin(qJ(2));
	t232 = sin(qJ(1));
	t246 = t231 * t232;
	t234 = cos(qJ(1));
	t245 = t231 * t234;
	t233 = cos(qJ(2));
	t244 = t232 * t233;
	t243 = t233 * t234;
	t228 = sin(pkin(6));
	t242 = qJD(1) * t228;
	t241 = qJD(2) * t231;
	t230 = cos(pkin(6));
	t240 = t230 * t246;
	t239 = t232 * t242;
	t238 = t234 * t242;
	t237 = t228 * t241;
	t236 = -t230 * t244 - t245;
	t235 = -t230 * t245 - t244;
	t229 = cos(pkin(11));
	t227 = sin(pkin(11));
	t224 = -qJD(1) * t240 - t232 * t241 + (qJD(2) * t230 + qJD(1)) * t243;
	t223 = t236 * qJD(1) + t235 * qJD(2);
	t222 = t235 * qJD(1) + t236 * qJD(2);
	t221 = (t240 - t243) * qJD(2) + (-t230 * t243 + t246) * qJD(1);
	t1 = [-t224 * t229 - t227 * t239, t221 * t229, 0, 0, 0, 0; t222 * t229 + t227 * t238, t223 * t229, 0, 0, 0, 0; 0, -t229 * t237, 0, 0, 0, 0; t224 * t227 - t229 * t239, -t221 * t227, 0, 0, 0, 0; -t222 * t227 + t229 * t238, -t223 * t227, 0, 0, 0, 0; 0, t227 * t237, 0, 0, 0, 0; t223, t222, 0, 0, 0, 0; -t221, t224, 0, 0, 0, 0; 0, t228 * qJD(2) * t233, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:19:00
	% EndTime: 2019-10-10 10:19:01
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (144->36), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->38)
	t303 = cos(pkin(6));
	t305 = sin(qJ(1));
	t304 = sin(qJ(2));
	t324 = t305 * t304;
	t315 = t303 * t324;
	t319 = qJD(2) * t304;
	t306 = cos(qJ(2));
	t307 = cos(qJ(1));
	t321 = t307 * t306;
	t291 = -qJD(1) * t315 - t305 * t319 + (qJD(2) * t303 + qJD(1)) * t321;
	t322 = t307 * t304;
	t323 = t305 * t306;
	t293 = t303 * t322 + t323;
	t301 = pkin(11) + qJ(4);
	t299 = sin(t301);
	t300 = cos(t301);
	t302 = sin(pkin(6));
	t320 = qJD(1) * t302;
	t314 = t305 * t320;
	t325 = t302 * t307;
	t328 = (-t293 * t300 + t299 * t325) * qJD(4) - t291 * t299 + t300 * t314;
	t327 = t302 * t304;
	t326 = t302 * t305;
	t318 = qJD(4) * t299;
	t317 = qJD(4) * t300;
	t316 = qJD(4) * t306;
	t313 = t307 * t320;
	t312 = t302 * qJD(2) * t306;
	t292 = t303 * t321 - t324;
	t294 = -t303 * t323 - t322;
	t310 = t315 - t321;
	t308 = -t291 * t300 + t317 * t325 + (qJD(4) * t293 - t314) * t299;
	t290 = qJD(1) * t294 - qJD(2) * t293;
	t289 = -qJD(1) * t293 + qJD(2) * t294;
	t288 = -qJD(1) * t292 + qJD(2) * t310;
	t287 = t299 * t313 + t289 * t300 + (t299 * t310 + t300 * t326) * qJD(4);
	t286 = t300 * t313 - t289 * t299 + (-t299 * t326 + t300 * t310) * qJD(4);
	t1 = [t308, t288 * t300 - t294 * t318, 0, t286, 0, 0; t287, t290 * t300 - t292 * t318, 0, t328, 0, 0; 0, (-t299 * t316 - t300 * t319) * t302, 0, -t299 * t312 + (-t299 * t303 - t300 * t327) * qJD(4), 0, 0; -t328, -t288 * t299 - t294 * t317, 0, -t287, 0, 0; t286, -t290 * t299 - t292 * t317, 0, t308, 0, 0; 0, (t299 * t319 - t300 * t316) * t302, 0, -t300 * t312 + (t299 * t327 - t300 * t303) * qJD(4), 0, 0; t290, t289, 0, 0, 0, 0; -t288, t291, 0, 0, 0, 0; 0, t312, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:19:02
	% EndTime: 2019-10-10 10:19:02
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (235->49), mult. (518->103), div. (0->0), fcn. (536->10), ass. (0->48)
	t412 = sin(qJ(1));
	t410 = cos(pkin(6));
	t419 = qJD(2) * t410 + qJD(1);
	t411 = sin(qJ(2));
	t435 = t412 * t411;
	t425 = t410 * t435;
	t430 = qJD(2) * t411;
	t413 = cos(qJ(2));
	t414 = cos(qJ(1));
	t432 = t414 * t413;
	t391 = -qJD(1) * t425 - t412 * t430 + t419 * t432;
	t433 = t414 * t411;
	t434 = t412 * t413;
	t394 = t410 * t433 + t434;
	t406 = pkin(11) + qJ(4);
	t404 = sin(t406);
	t405 = cos(t406);
	t408 = sin(pkin(6));
	t431 = qJD(1) * t408;
	t423 = t412 * t431;
	t436 = t408 * t414;
	t387 = (t394 * t404 + t405 * t436) * qJD(4) - t391 * t405 - t404 * t423;
	t439 = t405 * t411;
	t438 = t408 * t411;
	t437 = t408 * t412;
	t429 = qJD(2) * t413;
	t428 = qJD(4) * t404;
	t427 = qJD(4) * t405;
	t426 = qJD(4) * t413;
	t424 = t410 * t432;
	t422 = t414 * t431;
	t421 = t408 * t429;
	t420 = t404 * t426;
	t395 = -t410 * t434 - t433;
	t388 = -qJD(1) * t424 - t414 * t429 + t419 * t435;
	t417 = -t388 * t405 + t395 * t428;
	t390 = t395 * qJD(1) - t394 * qJD(2);
	t393 = t424 - t435;
	t416 = -t390 * t405 + t393 * t428;
	t386 = -t391 * t404 - t394 * t427 + t405 * t423 + t428 * t436;
	t409 = cos(pkin(12));
	t407 = sin(pkin(12));
	t396 = -t425 + t432;
	t392 = -t404 * t421 + (-t404 * t410 - t405 * t438) * qJD(4);
	t389 = -t394 * qJD(1) + t395 * qJD(2);
	t385 = t404 * t422 + t389 * t405 + (-t396 * t404 + t405 * t437) * qJD(4);
	t384 = t389 * t404 - t405 * t422 + (t396 * t405 + t404 * t437) * qJD(4);
	t1 = [t387 * t409 + t390 * t407, t389 * t407 - t417 * t409, 0, -t384 * t409, 0, 0; t385 * t409 - t388 * t407, t391 * t407 - t416 * t409, 0, t386 * t409, 0, 0; 0, (-t409 * t420 + (t407 * t413 - t409 * t439) * qJD(2)) * t408, 0, t392 * t409, 0, 0; -t387 * t407 + t390 * t409, t389 * t409 + t417 * t407, 0, t384 * t407, 0, 0; -t385 * t407 - t388 * t409, t391 * t409 + t416 * t407, 0, -t386 * t407, 0, 0; 0, (t407 * t420 + (t407 * t439 + t409 * t413) * qJD(2)) * t408, 0, -t392 * t407, 0, 0; t386, t388 * t404 + t395 * t427, 0, t385, 0, 0; t384, t390 * t404 + t393 * t427, 0, -t387, 0, 0; 0, (-t404 * t430 + t405 * t426) * t408, 0, t405 * t421 + (-t404 * t438 + t405 * t410) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:19:02
	% EndTime: 2019-10-10 10:19:03
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (506->75), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->68)
	t491 = sin(qJ(1));
	t489 = cos(pkin(6));
	t502 = qJD(2) * t489 + qJD(1);
	t490 = sin(qJ(2));
	t523 = t491 * t490;
	t510 = t489 * t523;
	t518 = qJD(2) * t490;
	t492 = cos(qJ(2));
	t493 = cos(qJ(1));
	t520 = t493 * t492;
	t460 = -qJD(1) * t510 - t491 * t518 + t502 * t520;
	t487 = pkin(11) + qJ(4);
	t483 = sin(t487);
	t485 = cos(t487);
	t521 = t493 * t490;
	t522 = t491 * t492;
	t471 = t489 * t521 + t522;
	t488 = sin(pkin(6));
	t524 = t488 * t493;
	t497 = t471 * t483 + t485 * t524;
	t519 = qJD(1) * t488;
	t507 = t491 * t519;
	t456 = t497 * qJD(4) - t460 * t485 - t483 * t507;
	t472 = t489 * t522 + t521;
	t459 = t472 * qJD(1) + t471 * qJD(2);
	t509 = t483 * t524;
	t465 = -t471 * t485 + t509;
	t508 = t489 * t520;
	t470 = -t508 + t523;
	t486 = pkin(12) + qJ(6);
	t482 = sin(t486);
	t484 = cos(t486);
	t536 = t456 * t484 - t459 * t482 + (-t465 * t482 - t470 * t484) * qJD(6);
	t535 = (t465 * t484 - t470 * t482) * qJD(6) + t456 * t482 + t459 * t484;
	t514 = qJD(4) * t492;
	t532 = (qJD(2) * t485 - qJD(6)) * t490 + t483 * t514;
	t527 = t488 * t490;
	t526 = t488 * t491;
	t525 = t488 * t492;
	t517 = qJD(2) * t492;
	t516 = qJD(4) * t483;
	t515 = qJD(4) * t485;
	t513 = qJD(6) * t482;
	t512 = qJD(6) * t484;
	t511 = qJD(6) * t485;
	t506 = t493 * t519;
	t505 = t488 * t518;
	t504 = t488 * t517;
	t458 = -t471 * qJD(1) - t472 * qJD(2);
	t500 = t472 * t511 + t458;
	t499 = t470 * t511 + t460;
	t498 = (qJD(2) - t511) * t492;
	t473 = -t510 + t520;
	t466 = -t473 * t483 + t485 * t526;
	t467 = t473 * t485 + t483 * t526;
	t469 = t489 * t483 + t485 * t527;
	t468 = -t483 * t527 + t489 * t485;
	t454 = qJD(4) * t509 - t460 * t483 - t471 * t515 + t485 * t507;
	t457 = -qJD(1) * t508 - t493 * t517 + t502 * t523;
	t495 = qJD(6) * t473 + t457 * t485 + t472 * t516;
	t494 = qJD(6) * t471 - t459 * t485 + t470 * t516;
	t462 = t468 * qJD(4) + t485 * t504;
	t461 = -t469 * qJD(4) - t483 * t504;
	t453 = t466 * qJD(4) + t458 * t485 + t483 * t506;
	t452 = t467 * qJD(4) + t458 * t483 - t485 * t506;
	t451 = t453 * t484 - t457 * t482 + (-t467 * t482 + t472 * t484) * qJD(6);
	t450 = -t453 * t482 - t457 * t484 + (-t467 * t484 - t472 * t482) * qJD(6);
	t1 = [t536, t482 * t500 + t484 * t495, 0, -t452 * t484 - t466 * t513, 0, t450; t451, t482 * t499 + t484 * t494, 0, t454 * t484 + t497 * t513, 0, t535; 0, (t482 * t498 - t532 * t484) * t488, 0, t461 * t484 - t468 * t513, 0, t484 * t505 - t462 * t482 + (-t469 * t484 + t482 * t525) * qJD(6); -t535, -t482 * t495 + t484 * t500, 0, t452 * t482 - t466 * t512, 0, -t451; t450, -t482 * t494 + t484 * t499, 0, -t454 * t482 + t497 * t512, 0, t536; 0, (t532 * t482 + t484 * t498) * t488, 0, -t461 * t482 - t468 * t512, 0, -t482 * t505 - t462 * t484 + (t469 * t482 + t484 * t525) * qJD(6); t454, t457 * t483 - t472 * t515, 0, t453, 0, 0; t452, -t459 * t483 - t470 * t515, 0, -t456, 0, 0; 0, (-t483 * t518 + t485 * t514) * t488, 0, t462, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end