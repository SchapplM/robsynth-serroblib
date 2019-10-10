% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:48
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRRP13_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP13_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:31
	% EndTime: 2019-10-10 10:48:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:31
	% EndTime: 2019-10-10 10:48:31
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
	% StartTime: 2019-10-10 10:48:31
	% EndTime: 2019-10-10 10:48:31
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
	% StartTime: 2019-10-10 10:48:32
	% EndTime: 2019-10-10 10:48:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->11), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t201 = sin(qJ(2));
	t202 = sin(qJ(1));
	t214 = t201 * t202;
	t204 = cos(qJ(1));
	t213 = t201 * t204;
	t203 = cos(qJ(2));
	t212 = t202 * t203;
	t211 = t203 * t204;
	t199 = sin(pkin(6));
	t210 = qJD(1) * t199;
	t209 = qJD(2) * t199;
	t200 = cos(pkin(6));
	t208 = t200 * t211 - t214;
	t207 = t200 * t212 + t213;
	t206 = t200 * t213 + t212;
	t205 = -t200 * t214 + t211;
	t198 = t205 * qJD(1) + t208 * qJD(2);
	t197 = t207 * qJD(1) + t206 * qJD(2);
	t196 = t206 * qJD(1) + t207 * qJD(2);
	t195 = t208 * qJD(1) + t205 * qJD(2);
	t1 = [-t202 * t210, 0, 0, 0, 0, 0; t204 * t210, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t198, t195, 0, 0, 0, 0; t196, t197, 0, 0, 0, 0; 0, t201 * t209, 0, 0, 0, 0; -t197, -t196, 0, 0, 0, 0; t195, t198, 0, 0, 0, 0; 0, t203 * t209, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:32
	% EndTime: 2019-10-10 10:48:32
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (95->37), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->36)
	t290 = cos(pkin(6));
	t292 = sin(qJ(2));
	t296 = cos(qJ(1));
	t310 = t296 * t292;
	t293 = sin(qJ(1));
	t295 = cos(qJ(2));
	t311 = t293 * t295;
	t283 = t290 * t310 + t311;
	t284 = t290 * t311 + t310;
	t280 = qJD(1) * t284 + qJD(2) * t283;
	t291 = sin(qJ(4));
	t294 = cos(qJ(4));
	t309 = t296 * t295;
	t304 = t290 * t309;
	t312 = t293 * t292;
	t299 = t304 - t312;
	t289 = sin(pkin(6));
	t308 = qJD(1) * t289;
	t303 = t293 * t308;
	t313 = t289 * t296;
	t316 = qJD(4) * (t291 * t299 + t294 * t313) + t280 * t294 - t291 * t303;
	t315 = t289 * t293;
	t314 = t289 * t295;
	t307 = qJD(2) * t295;
	t306 = qJD(4) * t291;
	t305 = qJD(4) * t294;
	t302 = t296 * t308;
	t301 = t289 * qJD(2) * t292;
	t285 = -t290 * t312 + t309;
	t297 = -t294 * t303 - t280 * t291 + (-t291 * t313 + t294 * t299) * qJD(4);
	t281 = qJD(1) * t285 + qJD(2) * t299;
	t279 = -qJD(1) * t283 - qJD(2) * t284;
	t278 = -qJD(1) * t304 - t296 * t307 + (qJD(2) * t290 + qJD(1)) * t312;
	t277 = t294 * t302 - t278 * t291 + (t284 * t294 - t291 * t315) * qJD(4);
	t276 = -t291 * t302 - t278 * t294 + (-t284 * t291 - t294 * t315) * qJD(4);
	t1 = [t297, t279 * t291 + t285 * t305, 0, t276, 0, 0; t277, t281 * t291 + t283 * t305, 0, t316, 0, 0; 0, (t291 * t307 + t292 * t305) * t289, 0, t294 * t301 + (-t290 * t294 + t291 * t314) * qJD(4), 0, 0; -t316, t279 * t294 - t285 * t306, 0, -t277, 0, 0; t276, t281 * t294 - t283 * t306, 0, t297, 0, 0; 0, (-t292 * t306 + t294 * t307) * t289, 0, -t291 * t301 + (t290 * t291 + t294 * t314) * qJD(4), 0, 0; -t281, t278, 0, 0, 0, 0; t279, -t280, 0, 0, 0, 0; 0, -t301, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:34
	% EndTime: 2019-10-10 10:48:35
	% DurationCPUTime: 0.62s
	% Computational Cost: add. (289->72), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->63)
	t459 = cos(pkin(6));
	t462 = sin(qJ(2));
	t467 = cos(qJ(1));
	t492 = t467 * t462;
	t463 = sin(qJ(1));
	t466 = cos(qJ(2));
	t493 = t463 * t466;
	t451 = t459 * t492 + t493;
	t452 = t459 * t493 + t492;
	t439 = t452 * qJD(1) + t451 * qJD(2);
	t491 = t467 * t466;
	t494 = t463 * t462;
	t450 = -t459 * t491 + t494;
	t461 = sin(qJ(4));
	t465 = cos(qJ(4));
	t458 = sin(pkin(6));
	t496 = t458 * t467;
	t445 = t450 * t465 + t461 * t496;
	t490 = qJD(1) * t458;
	t480 = t463 * t490;
	t432 = t445 * qJD(4) + t439 * t461 + t465 * t480;
	t482 = t459 * t494;
	t489 = qJD(2) * t462;
	t440 = -qJD(1) * t482 - t463 * t489 + (qJD(2) * t459 + qJD(1)) * t491;
	t481 = t465 * t496;
	t446 = -t450 * t461 + t481;
	t460 = sin(qJ(5));
	t464 = cos(qJ(5));
	t507 = -t432 * t464 + (-t446 * t460 - t451 * t464) * qJD(5) - t440 * t460;
	t506 = (t446 * t464 - t451 * t460) * qJD(5) - t432 * t460 + t440 * t464;
	t503 = (qJD(2) * t461 + qJD(5)) * t466;
	t498 = t458 * t463;
	t497 = t458 * t466;
	t495 = t462 * t464;
	t488 = qJD(2) * t466;
	t487 = qJD(4) * t461;
	t486 = qJD(4) * t465;
	t485 = qJD(5) * t460;
	t484 = qJD(5) * t461;
	t483 = qJD(5) * t464;
	t479 = t467 * t490;
	t478 = t458 * t489;
	t477 = t458 * t488;
	t476 = -qJD(2) - t484;
	t471 = t482 - t491;
	t437 = t450 * qJD(1) + t471 * qJD(2);
	t474 = t471 * t484 + t437;
	t473 = -t451 * t484 - t439;
	t444 = t452 * t461 + t465 * t498;
	t443 = t452 * t465 - t461 * t498;
	t448 = -t459 * t461 - t465 * t497;
	t472 = -t459 * t465 + t461 * t497;
	t438 = -t451 * qJD(1) - t452 * qJD(2);
	t469 = -qJD(5) * t452 + t438 * t461 - t471 * t486;
	t468 = -qJD(5) * t450 + t440 * t461 + t451 * t486;
	t431 = t439 * t465 + qJD(4) * t481 + (-qJD(4) * t450 - t480) * t461;
	t442 = t472 * qJD(4) + t465 * t478;
	t441 = t448 * qJD(4) + t461 * t478;
	t435 = t443 * qJD(4) - t437 * t461 + t465 * t479;
	t434 = t444 * qJD(4) + t437 * t465 + t461 * t479;
	t430 = t435 * t464 + t438 * t460 + (-t444 * t460 - t464 * t471) * qJD(5);
	t429 = -t435 * t460 + t438 * t464 + (-t444 * t464 + t460 * t471) * qJD(5);
	t1 = [t507, t474 * t460 + t469 * t464, 0, -t434 * t464 - t443 * t485, t429, 0; t430, t473 * t460 + t468 * t464, 0, t431 * t464 - t445 * t485, t506, 0; 0, (t464 * t503 + (t476 * t460 + t464 * t486) * t462) * t458, 0, t442 * t464 - t448 * t485, t464 * t477 - t441 * t460 + (-t458 * t460 * t462 + t464 * t472) * qJD(5), 0; -t506, -t469 * t460 + t474 * t464, 0, t434 * t460 - t443 * t483, -t430, 0; t429, -t468 * t460 + t473 * t464, 0, -t431 * t460 - t445 * t483, t507, 0; 0, (t476 * t495 + (-t462 * t486 - t503) * t460) * t458, 0, -t442 * t460 - t448 * t483, -t460 * t477 - t441 * t464 + (-t458 * t495 - t460 * t472) * qJD(5), 0; t431, -t438 * t465 - t471 * t487, 0, t435, 0, 0; t434, -t440 * t465 + t451 * t487, 0, t432, 0, 0; 0, (t462 * t487 - t465 * t488) * t458, 0, t441, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:48:35
	% EndTime: 2019-10-10 10:48:35
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (289->72), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->63)
	t481 = cos(pkin(6));
	t484 = sin(qJ(2));
	t489 = cos(qJ(1));
	t514 = t489 * t484;
	t485 = sin(qJ(1));
	t488 = cos(qJ(2));
	t515 = t485 * t488;
	t473 = t481 * t514 + t515;
	t474 = t481 * t515 + t514;
	t461 = t474 * qJD(1) + t473 * qJD(2);
	t513 = t489 * t488;
	t516 = t485 * t484;
	t472 = -t481 * t513 + t516;
	t483 = sin(qJ(4));
	t487 = cos(qJ(4));
	t480 = sin(pkin(6));
	t518 = t480 * t489;
	t467 = t472 * t487 + t483 * t518;
	t512 = qJD(1) * t480;
	t502 = t485 * t512;
	t454 = t467 * qJD(4) + t461 * t483 + t487 * t502;
	t504 = t481 * t516;
	t511 = qJD(2) * t484;
	t462 = -qJD(1) * t504 - t485 * t511 + (qJD(2) * t481 + qJD(1)) * t513;
	t503 = t487 * t518;
	t468 = -t472 * t483 + t503;
	t482 = sin(qJ(5));
	t486 = cos(qJ(5));
	t529 = -t454 * t486 + (-t468 * t482 - t473 * t486) * qJD(5) - t462 * t482;
	t528 = (t468 * t486 - t473 * t482) * qJD(5) - t454 * t482 + t462 * t486;
	t525 = (qJD(2) * t483 + qJD(5)) * t488;
	t520 = t480 * t485;
	t519 = t480 * t488;
	t517 = t484 * t486;
	t510 = qJD(2) * t488;
	t509 = qJD(4) * t483;
	t508 = qJD(4) * t487;
	t507 = qJD(5) * t482;
	t506 = qJD(5) * t483;
	t505 = qJD(5) * t486;
	t501 = t489 * t512;
	t500 = t480 * t511;
	t499 = t480 * t510;
	t498 = -qJD(2) - t506;
	t493 = t504 - t513;
	t459 = t472 * qJD(1) + t493 * qJD(2);
	t496 = t493 * t506 + t459;
	t495 = -t473 * t506 - t461;
	t466 = t474 * t483 + t487 * t520;
	t465 = t474 * t487 - t483 * t520;
	t470 = -t481 * t483 - t487 * t519;
	t494 = -t481 * t487 + t483 * t519;
	t460 = -t473 * qJD(1) - t474 * qJD(2);
	t491 = -qJD(5) * t474 + t460 * t483 - t493 * t508;
	t490 = -qJD(5) * t472 + t462 * t483 + t473 * t508;
	t453 = t461 * t487 + qJD(4) * t503 + (-qJD(4) * t472 - t502) * t483;
	t464 = t494 * qJD(4) + t487 * t500;
	t463 = t470 * qJD(4) + t483 * t500;
	t457 = t465 * qJD(4) - t459 * t483 + t487 * t501;
	t456 = t466 * qJD(4) + t459 * t487 + t483 * t501;
	t452 = t457 * t486 + t460 * t482 + (-t466 * t482 - t486 * t493) * qJD(5);
	t451 = -t457 * t482 + t460 * t486 + (-t466 * t486 + t482 * t493) * qJD(5);
	t1 = [t529, t496 * t482 + t491 * t486, 0, -t456 * t486 - t465 * t507, t451, 0; t452, t495 * t482 + t490 * t486, 0, t453 * t486 - t467 * t507, t528, 0; 0, (t486 * t525 + (t498 * t482 + t486 * t508) * t484) * t480, 0, t464 * t486 - t470 * t507, t486 * t499 - t463 * t482 + (-t480 * t482 * t484 + t486 * t494) * qJD(5), 0; -t528, -t491 * t482 + t496 * t486, 0, t456 * t482 - t465 * t505, -t452, 0; t451, -t490 * t482 + t495 * t486, 0, -t453 * t482 - t467 * t505, t529, 0; 0, (t498 * t517 + (-t484 * t508 - t525) * t482) * t480, 0, -t464 * t482 - t470 * t505, -t482 * t499 - t463 * t486 + (-t480 * t517 - t482 * t494) * qJD(5), 0; t453, -t460 * t487 - t493 * t509, 0, t457, 0, 0; t456, -t462 * t487 + t473 * t509, 0, t454, 0, 0; 0, (t484 * t509 - t487 * t510) * t480, 0, t463, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end