% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:43
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:07
	% EndTime: 2019-10-10 09:43:07
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
	% StartTime: 2019-10-10 09:43:07
	% EndTime: 2019-10-10 09:43:07
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (26->12), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t205 = sin(qJ(2));
	t206 = sin(qJ(1));
	t218 = t205 * t206;
	t208 = cos(qJ(1));
	t217 = t205 * t208;
	t207 = cos(qJ(2));
	t216 = t206 * t207;
	t215 = t207 * t208;
	t203 = sin(pkin(6));
	t214 = qJD(1) * t203;
	t213 = qJD(2) * t203;
	t204 = cos(pkin(6));
	t212 = t204 * t215 - t218;
	t211 = -t204 * t216 - t217;
	t210 = -t204 * t217 - t216;
	t209 = t204 * t218 - t215;
	t202 = -t209 * qJD(1) + t212 * qJD(2);
	t201 = t211 * qJD(1) + t210 * qJD(2);
	t200 = t210 * qJD(1) + t211 * qJD(2);
	t199 = -t212 * qJD(1) + t209 * qJD(2);
	t1 = [-t202, t199, 0, 0, 0, 0; t200, t201, 0, 0, 0, 0; 0, -t205 * t213, 0, 0, 0, 0; -t206 * t214, 0, 0, 0, 0, 0; t208 * t214, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t201, t200, 0, 0, 0, 0; -t199, t202, 0, 0, 0, 0; 0, t207 * t213, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:07
	% EndTime: 2019-10-10 09:43:07
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->12), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t175 = sin(qJ(2));
	t176 = sin(qJ(1));
	t188 = t175 * t176;
	t178 = cos(qJ(1));
	t187 = t175 * t178;
	t177 = cos(qJ(2));
	t186 = t176 * t177;
	t185 = t177 * t178;
	t173 = sin(pkin(6));
	t184 = qJD(1) * t173;
	t183 = qJD(2) * t173;
	t174 = cos(pkin(6));
	t182 = t174 * t185 - t188;
	t181 = -t174 * t186 - t187;
	t180 = -t174 * t187 - t186;
	t179 = t174 * t188 - t185;
	t172 = -t179 * qJD(1) + t182 * qJD(2);
	t171 = t181 * qJD(1) + t180 * qJD(2);
	t170 = t180 * qJD(1) + t181 * qJD(2);
	t169 = -t182 * qJD(1) + t179 * qJD(2);
	t1 = [-t172, t169, 0, 0, 0, 0; t170, t171, 0, 0, 0, 0; 0, -t175 * t183, 0, 0, 0, 0; t171, t170, 0, 0, 0, 0; -t169, t172, 0, 0, 0, 0; 0, t177 * t183, 0, 0, 0, 0; t176 * t184, 0, 0, 0, 0, 0; -t178 * t184, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:07
	% EndTime: 2019-10-10 09:43:08
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (95->36), mult. (310->71), div. (0->0), fcn. (322->8), ass. (0->37)
	t290 = cos(pkin(6));
	t293 = sin(qJ(1));
	t292 = sin(qJ(2));
	t314 = t293 * t292;
	t305 = t290 * t314;
	t309 = qJD(2) * t292;
	t295 = cos(qJ(2));
	t296 = cos(qJ(1));
	t311 = t296 * t295;
	t281 = -qJD(1) * t305 - t293 * t309 + (qJD(2) * t290 + qJD(1)) * t311;
	t312 = t296 * t292;
	t313 = t293 * t295;
	t283 = t290 * t312 + t313;
	t291 = sin(qJ(5));
	t294 = cos(qJ(5));
	t289 = sin(pkin(6));
	t310 = qJD(1) * t289;
	t304 = t293 * t310;
	t315 = t289 * t296;
	t318 = (t283 * t294 + t291 * t315) * qJD(5) + t281 * t291 + t294 * t304;
	t317 = t289 * t291;
	t316 = t289 * t294;
	t308 = qJD(5) * t291;
	t307 = qJD(5) * t294;
	t306 = qJD(5) * t295;
	t303 = t296 * t310;
	t302 = t289 * qJD(2) * t295;
	t282 = t290 * t311 - t314;
	t300 = t290 * t313 + t312;
	t299 = t305 - t311;
	t297 = -t281 * t294 + t291 * t304 + (t283 * t291 - t294 * t315) * qJD(5);
	t280 = t300 * qJD(1) + t283 * qJD(2);
	t279 = t283 * qJD(1) + t300 * qJD(2);
	t278 = -t282 * qJD(1) + t299 * qJD(2);
	t277 = -t291 * t303 - t279 * t294 + (t291 * t299 - t293 * t316) * qJD(5);
	t276 = -t294 * t303 + t279 * t291 + (t293 * t317 + t294 * t299) * qJD(5);
	t1 = [t297, t278 * t294 + t300 * t308, 0, 0, t276, 0; t277, -t280 * t294 - t282 * t308, 0, 0, -t318, 0; 0, (-t291 * t306 - t294 * t309) * t289, 0, 0, -t291 * t302 + (t290 * t291 - t292 * t316) * qJD(5), 0; t318, -t278 * t291 + t300 * t307, 0, 0, -t277, 0; t276, t280 * t291 - t282 * t307, 0, 0, t297, 0; 0, (t291 * t309 - t294 * t306) * t289, 0, 0, -t294 * t302 + (t290 * t294 + t292 * t317) * qJD(5), 0; t280, t279, 0, 0, 0, 0; t278, -t281, 0, 0, 0, 0; 0, -t302, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:09
	% EndTime: 2019-10-10 09:43:10
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (289->71), mult. (898->138), div. (0->0), fcn. (964->10), ass. (0->63)
	t448 = cos(pkin(6));
	t452 = sin(qJ(1));
	t451 = sin(qJ(2));
	t484 = t452 * t451;
	t471 = t448 * t484;
	t479 = qJD(2) * t451;
	t455 = cos(qJ(2));
	t456 = cos(qJ(1));
	t481 = t456 * t455;
	t428 = -qJD(1) * t471 - t452 * t479 + (qJD(2) * t448 + qJD(1)) * t481;
	t482 = t456 * t451;
	t483 = t452 * t455;
	t439 = t448 * t482 + t483;
	t450 = sin(qJ(5));
	t454 = cos(qJ(5));
	t447 = sin(pkin(6));
	t485 = t447 * t456;
	t431 = -t439 * t450 + t454 * t485;
	t480 = qJD(1) * t447;
	t470 = t452 * t480;
	t423 = t431 * qJD(5) + t428 * t454 - t450 * t470;
	t461 = t448 * t483 + t482;
	t427 = t461 * qJD(1) + t439 * qJD(2);
	t432 = t439 * t454 + t450 * t485;
	t438 = -t448 * t481 + t484;
	t449 = sin(qJ(6));
	t453 = cos(qJ(6));
	t497 = -t423 * t453 + t427 * t449 + (t432 * t449 + t438 * t453) * qJD(6);
	t496 = (t432 * t453 - t438 * t449) * qJD(6) + t423 * t449 + t427 * t453;
	t476 = qJD(5) * t455;
	t493 = (qJD(2) * t454 + qJD(6)) * t451 + t450 * t476;
	t488 = t447 * t450;
	t487 = t447 * t454;
	t486 = t447 * t455;
	t478 = qJD(5) * t450;
	t477 = qJD(5) * t454;
	t475 = qJD(6) * t449;
	t474 = qJD(6) * t453;
	t473 = qJD(6) * t454;
	t472 = t452 * t488;
	t469 = t456 * t480;
	t468 = t447 * t479;
	t467 = qJD(2) * t486;
	t426 = t439 * qJD(1) + t461 * qJD(2);
	t464 = t461 * t473 + t426;
	t463 = t438 * t473 - t428;
	t462 = (-qJD(2) - t473) * t455;
	t460 = t471 - t481;
	t434 = t450 * t460 - t452 * t487;
	t437 = -t448 * t450 + t451 * t487;
	t436 = -t448 * t454 - t451 * t488;
	t425 = t438 * qJD(1) + t460 * qJD(2);
	t458 = -qJD(6) * t460 - t425 * t454 - t461 * t478;
	t457 = -qJD(6) * t439 - t427 * t454 + t438 * t478;
	t422 = -t432 * qJD(5) - t428 * t450 - t454 * t470;
	t435 = -t454 * t460 - t472;
	t430 = t436 * qJD(5) + t454 * t467;
	t429 = -t437 * qJD(5) - t450 * t467;
	t421 = t434 * qJD(5) - t426 * t454 - t450 * t469;
	t420 = -t426 * t450 - qJD(5) * t472 + (-qJD(5) * t460 + t469) * t454;
	t419 = t421 * t453 + t425 * t449 + (-t435 * t449 - t453 * t461) * qJD(6);
	t418 = -t421 * t449 + t425 * t453 + (-t435 * t453 + t449 * t461) * qJD(6);
	t1 = [t497, t464 * t449 - t458 * t453, 0, 0, -t420 * t453 - t434 * t475, t418; t419, t463 * t449 + t457 * t453, 0, 0, t422 * t453 - t431 * t475, -t496; 0, (t449 * t462 - t493 * t453) * t447, 0, 0, t429 * t453 - t436 * t475, -t453 * t468 - t430 * t449 + (-t437 * t453 - t449 * t486) * qJD(6); t496, t458 * t449 + t464 * t453, 0, 0, t420 * t449 - t434 * t474, -t419; t418, -t457 * t449 + t463 * t453, 0, 0, -t422 * t449 - t431 * t474, t497; 0, (t493 * t449 + t453 * t462) * t447, 0, 0, -t429 * t449 - t436 * t474, t449 * t468 - t430 * t453 + (t437 * t449 - t453 * t486) * qJD(6); t422, t425 * t450 - t461 * t477, 0, 0, t421, 0; t420, -t427 * t450 - t438 * t477, 0, 0, t423, 0; 0, (-t450 * t479 + t454 * t476) * t447, 0, 0, t430, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end