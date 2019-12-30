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
% Datum: 2019-12-29 19:30
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
	% StartTime: 2019-12-29 19:30:31
	% EndTime: 2019-12-29 19:30:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:30
	% EndTime: 2019-12-29 19:30:30
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:36
	% EndTime: 2019-12-29 19:30:36
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-12-29 19:30:31
	% EndTime: 2019-12-29 19:30:31
	% DurationCPUTime: 0.16s
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
	% StartTime: 2019-12-29 19:30:32
	% EndTime: 2019-12-29 19:30:32
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (95->37), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->36)
	t289 = cos(pkin(5));
	t291 = sin(qJ(2));
	t295 = cos(qJ(1));
	t309 = t295 * t291;
	t292 = sin(qJ(1));
	t294 = cos(qJ(2));
	t310 = t292 * t294;
	t282 = t289 * t309 + t310;
	t283 = t289 * t310 + t309;
	t279 = t283 * qJD(1) + t282 * qJD(2);
	t290 = sin(qJ(4));
	t293 = cos(qJ(4));
	t308 = t295 * t294;
	t303 = t289 * t308;
	t311 = t292 * t291;
	t298 = t303 - t311;
	t288 = sin(pkin(5));
	t307 = qJD(1) * t288;
	t302 = t292 * t307;
	t312 = t288 * t295;
	t315 = (t290 * t298 + t293 * t312) * qJD(4) + t279 * t293 - t290 * t302;
	t314 = t288 * t292;
	t313 = t288 * t294;
	t306 = qJD(2) * t294;
	t305 = qJD(4) * t290;
	t304 = qJD(4) * t293;
	t301 = t295 * t307;
	t300 = t288 * qJD(2) * t291;
	t284 = -t289 * t311 + t308;
	t296 = -t293 * t302 - t279 * t290 + (-t290 * t312 + t293 * t298) * qJD(4);
	t280 = t284 * qJD(1) + t298 * qJD(2);
	t278 = -t282 * qJD(1) - t283 * qJD(2);
	t277 = -qJD(1) * t303 - t295 * t306 + (qJD(2) * t289 + qJD(1)) * t311;
	t276 = t293 * t301 - t277 * t290 + (t283 * t293 - t290 * t314) * qJD(4);
	t275 = -t290 * t301 - t277 * t293 + (-t283 * t290 - t293 * t314) * qJD(4);
	t1 = [t296, t278 * t290 + t284 * t304, 0, t275, 0; t276, t280 * t290 + t282 * t304, 0, t315, 0; 0, (t290 * t306 + t291 * t304) * t288, 0, t293 * t300 + (-t289 * t293 + t290 * t313) * qJD(4), 0; -t315, t278 * t293 - t284 * t305, 0, -t276, 0; t275, t280 * t293 - t282 * t305, 0, t296, 0; 0, (-t291 * t305 + t293 * t306) * t288, 0, -t290 * t300 + (t289 * t290 + t293 * t313) * qJD(4), 0; -t280, t277, 0, 0, 0; t278, -t279, 0, 0, 0; 0, -t300, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:35
	% EndTime: 2019-12-29 19:30:36
	% DurationCPUTime: 0.88s
	% Computational Cost: add. (289->72), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->63)
	t461 = cos(pkin(5));
	t464 = sin(qJ(2));
	t469 = cos(qJ(1));
	t494 = t469 * t464;
	t465 = sin(qJ(1));
	t468 = cos(qJ(2));
	t495 = t465 * t468;
	t453 = t461 * t494 + t495;
	t454 = t461 * t495 + t494;
	t441 = t454 * qJD(1) + t453 * qJD(2);
	t493 = t469 * t468;
	t496 = t465 * t464;
	t452 = -t461 * t493 + t496;
	t463 = sin(qJ(4));
	t467 = cos(qJ(4));
	t460 = sin(pkin(5));
	t498 = t460 * t469;
	t447 = t452 * t467 + t463 * t498;
	t492 = qJD(1) * t460;
	t482 = t465 * t492;
	t434 = t447 * qJD(4) + t441 * t463 + t467 * t482;
	t484 = t461 * t496;
	t491 = qJD(2) * t464;
	t442 = -qJD(1) * t484 - t465 * t491 + (qJD(2) * t461 + qJD(1)) * t493;
	t483 = t467 * t498;
	t448 = -t452 * t463 + t483;
	t462 = sin(qJ(5));
	t466 = cos(qJ(5));
	t509 = -t434 * t466 + (-t448 * t462 - t453 * t466) * qJD(5) - t442 * t462;
	t508 = (t448 * t466 - t453 * t462) * qJD(5) - t434 * t462 + t442 * t466;
	t505 = (qJD(2) * t463 + qJD(5)) * t468;
	t500 = t460 * t465;
	t499 = t460 * t468;
	t497 = t464 * t466;
	t490 = qJD(2) * t468;
	t489 = qJD(4) * t463;
	t488 = qJD(4) * t467;
	t487 = qJD(5) * t462;
	t486 = qJD(5) * t463;
	t485 = qJD(5) * t466;
	t481 = t469 * t492;
	t480 = t460 * t491;
	t479 = t460 * t490;
	t478 = -qJD(2) - t486;
	t473 = t484 - t493;
	t439 = t452 * qJD(1) + t473 * qJD(2);
	t476 = t473 * t486 + t439;
	t475 = -t453 * t486 - t441;
	t446 = t454 * t463 + t467 * t500;
	t445 = t454 * t467 - t463 * t500;
	t450 = -t461 * t463 - t467 * t499;
	t474 = -t461 * t467 + t463 * t499;
	t440 = -t453 * qJD(1) - t454 * qJD(2);
	t471 = -qJD(5) * t454 + t440 * t463 - t473 * t488;
	t470 = -qJD(5) * t452 + t442 * t463 + t453 * t488;
	t433 = t441 * t467 + qJD(4) * t483 + (-qJD(4) * t452 - t482) * t463;
	t444 = t474 * qJD(4) + t467 * t480;
	t443 = t450 * qJD(4) + t463 * t480;
	t437 = t445 * qJD(4) - t439 * t463 + t467 * t481;
	t436 = t446 * qJD(4) + t439 * t467 + t463 * t481;
	t432 = t437 * t466 + t440 * t462 + (-t446 * t462 - t466 * t473) * qJD(5);
	t431 = -t437 * t462 + t440 * t466 + (-t446 * t466 + t462 * t473) * qJD(5);
	t1 = [t509, t476 * t462 + t471 * t466, 0, -t436 * t466 - t445 * t487, t431; t432, t475 * t462 + t470 * t466, 0, t433 * t466 - t447 * t487, t508; 0, (t466 * t505 + (t478 * t462 + t466 * t488) * t464) * t460, 0, t444 * t466 - t450 * t487, t466 * t479 - t443 * t462 + (-t460 * t462 * t464 + t466 * t474) * qJD(5); -t508, -t471 * t462 + t476 * t466, 0, t436 * t462 - t445 * t485, -t432; t431, -t470 * t462 + t475 * t466, 0, -t433 * t462 - t447 * t485, t509; 0, (t478 * t497 + (-t464 * t488 - t505) * t462) * t460, 0, -t444 * t462 - t450 * t485, -t462 * t479 - t443 * t466 + (-t460 * t497 - t462 * t474) * qJD(5); t433, -t440 * t467 - t473 * t489, 0, t437, 0; t436, -t442 * t467 + t453 * t489, 0, t434, 0; 0, (t464 * t489 - t467 * t490) * t460, 0, t443, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end