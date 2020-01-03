% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRR16_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR16_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:28
	% EndTime: 2019-12-31 20:48:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:28
	% EndTime: 2019-12-31 20:48:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:28
	% EndTime: 2019-12-31 20:48:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(5));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(5));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0; -t153, -t154, 0, 0, 0; 0, -t158 * t166, 0, 0, 0; t154, t153, 0, 0, 0; t152, t155, 0, 0, 0; 0, -t160 * t166, 0, 0, 0; -t159 * t167, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:28
	% EndTime: 2019-12-31 20:48:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (25->11), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t201 = sin(qJ(2));
	t202 = sin(qJ(1));
	t214 = t201 * t202;
	t204 = cos(qJ(1));
	t213 = t201 * t204;
	t203 = cos(qJ(2));
	t212 = t202 * t203;
	t211 = t203 * t204;
	t199 = sin(pkin(5));
	t210 = qJD(1) * t199;
	t209 = qJD(2) * t199;
	t200 = cos(pkin(5));
	t208 = t200 * t211 - t214;
	t207 = t200 * t212 + t213;
	t206 = t200 * t213 + t212;
	t205 = -t200 * t214 + t211;
	t198 = t205 * qJD(1) + t208 * qJD(2);
	t197 = t207 * qJD(1) + t206 * qJD(2);
	t196 = t206 * qJD(1) + t207 * qJD(2);
	t195 = t208 * qJD(1) + t205 * qJD(2);
	t1 = [-t202 * t210, 0, 0, 0, 0; t204 * t210, 0, 0, 0, 0; 0, 0, 0, 0, 0; t198, t195, 0, 0, 0; t196, t197, 0, 0, 0; 0, t201 * t209, 0, 0, 0; -t197, -t196, 0, 0, 0; t195, t198, 0, 0, 0; 0, t203 * t209, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:29
	% EndTime: 2019-12-31 20:48:29
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (95->37), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->36)
	t290 = cos(pkin(5));
	t292 = sin(qJ(2));
	t296 = cos(qJ(1));
	t310 = t296 * t292;
	t293 = sin(qJ(1));
	t295 = cos(qJ(2));
	t311 = t293 * t295;
	t283 = t290 * t310 + t311;
	t284 = t290 * t311 + t310;
	t280 = t284 * qJD(1) + t283 * qJD(2);
	t291 = sin(qJ(4));
	t294 = cos(qJ(4));
	t309 = t296 * t295;
	t304 = t290 * t309;
	t312 = t293 * t292;
	t299 = t304 - t312;
	t289 = sin(pkin(5));
	t308 = qJD(1) * t289;
	t303 = t293 * t308;
	t313 = t289 * t296;
	t316 = (t291 * t299 + t294 * t313) * qJD(4) + t280 * t294 - t291 * t303;
	t315 = t289 * t293;
	t314 = t289 * t295;
	t307 = qJD(2) * t295;
	t306 = qJD(4) * t291;
	t305 = qJD(4) * t294;
	t302 = t296 * t308;
	t301 = t289 * qJD(2) * t292;
	t285 = -t290 * t312 + t309;
	t297 = -t294 * t303 - t280 * t291 + (-t291 * t313 + t294 * t299) * qJD(4);
	t281 = t285 * qJD(1) + t299 * qJD(2);
	t279 = -t283 * qJD(1) - t284 * qJD(2);
	t278 = -qJD(1) * t304 - t296 * t307 + (qJD(2) * t290 + qJD(1)) * t312;
	t277 = t294 * t302 - t278 * t291 + (t284 * t294 - t291 * t315) * qJD(4);
	t276 = -t291 * t302 - t278 * t294 + (-t284 * t291 - t294 * t315) * qJD(4);
	t1 = [t297, t279 * t291 + t285 * t305, 0, t276, 0; t277, t281 * t291 + t283 * t305, 0, t316, 0; 0, (t291 * t307 + t292 * t305) * t289, 0, t294 * t301 + (-t290 * t294 + t291 * t314) * qJD(4), 0; -t316, t279 * t294 - t285 * t306, 0, -t277, 0; t276, t281 * t294 - t283 * t306, 0, t297, 0; 0, (-t292 * t306 + t294 * t307) * t289, 0, -t291 * t301 + (t290 * t291 + t294 * t314) * qJD(4), 0; -t281, t278, 0, 0, 0; t279, -t280, 0, 0, 0; 0, -t301, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:48:30
	% EndTime: 2019-12-31 20:48:31
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (289->72), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->63)
	t459 = cos(pkin(5));
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
	t458 = sin(pkin(5));
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
	t1 = [t507, t474 * t460 + t469 * t464, 0, -t434 * t464 - t443 * t485, t429; t430, t473 * t460 + t468 * t464, 0, t431 * t464 - t445 * t485, t506; 0, (t464 * t503 + (t476 * t460 + t464 * t486) * t462) * t458, 0, t442 * t464 - t448 * t485, t464 * t477 - t441 * t460 + (-t458 * t460 * t462 + t464 * t472) * qJD(5); -t506, -t469 * t460 + t474 * t464, 0, t434 * t460 - t443 * t483, -t430; t429, -t468 * t460 + t473 * t464, 0, -t431 * t460 - t445 * t483, t507; 0, (t476 * t495 + (-t462 * t486 - t503) * t460) * t458, 0, -t442 * t460 - t448 * t483, -t460 * t477 - t441 * t464 + (-t458 * t495 - t460 * t472) * qJD(5); t431, -t438 * t465 - t471 * t487, 0, t435, 0; t434, -t440 * t465 + t451 * t487, 0, t432, 0; 0, (t462 * t487 - t465 * t488) * t458, 0, t441, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end