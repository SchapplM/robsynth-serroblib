% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR7
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
% Datum: 2019-10-10 09:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
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
	% StartTime: 2019-10-10 09:46:41
	% EndTime: 2019-10-10 09:46:41
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
	% StartTime: 2019-10-10 09:46:41
	% EndTime: 2019-10-10 09:46:41
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
	% StartTime: 2019-10-10 09:46:41
	% EndTime: 2019-10-10 09:46:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->11), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t176 = sin(qJ(2));
	t177 = sin(qJ(1));
	t189 = t176 * t177;
	t179 = cos(qJ(1));
	t188 = t176 * t179;
	t178 = cos(qJ(2));
	t187 = t177 * t178;
	t186 = t178 * t179;
	t174 = sin(pkin(6));
	t185 = qJD(1) * t174;
	t184 = qJD(2) * t174;
	t175 = cos(pkin(6));
	t183 = t175 * t186 - t189;
	t182 = t175 * t187 + t188;
	t181 = t175 * t188 + t187;
	t180 = -t175 * t189 + t186;
	t173 = t180 * qJD(1) + t183 * qJD(2);
	t172 = t182 * qJD(1) + t181 * qJD(2);
	t171 = t181 * qJD(1) + t182 * qJD(2);
	t170 = t183 * qJD(1) + t180 * qJD(2);
	t1 = [-t172, -t171, 0, 0, 0, 0; t170, t173, 0, 0, 0, 0; 0, t178 * t184, 0, 0, 0, 0; t173, t170, 0, 0, 0, 0; t171, t172, 0, 0, 0, 0; 0, t176 * t184, 0, 0, 0, 0; t177 * t185, 0, 0, 0, 0, 0; -t179 * t185, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:42
	% EndTime: 2019-10-10 09:46:42
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (95->37), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->36)
	t288 = cos(pkin(6));
	t290 = sin(qJ(2));
	t294 = cos(qJ(1));
	t308 = t294 * t290;
	t291 = sin(qJ(1));
	t293 = cos(qJ(2));
	t309 = t291 * t293;
	t280 = t288 * t308 + t309;
	t281 = t288 * t309 + t308;
	t277 = t281 * qJD(1) + t280 * qJD(2);
	t289 = sin(qJ(5));
	t292 = cos(qJ(5));
	t307 = t294 * t293;
	t302 = t288 * t307;
	t310 = t291 * t290;
	t297 = t302 - t310;
	t287 = sin(pkin(6));
	t306 = qJD(1) * t287;
	t301 = t291 * t306;
	t311 = t287 * t294;
	t314 = (t289 * t311 - t292 * t297) * qJD(5) + t277 * t289 + t292 * t301;
	t313 = t287 * t291;
	t312 = t287 * t293;
	t305 = qJD(2) * t293;
	t304 = qJD(5) * t289;
	t303 = qJD(5) * t292;
	t300 = t294 * t306;
	t299 = t287 * qJD(2) * t290;
	t282 = -t288 * t310 + t307;
	t295 = -t277 * t292 + t289 * t301 + (-t289 * t297 - t292 * t311) * qJD(5);
	t278 = t282 * qJD(1) + t297 * qJD(2);
	t276 = -t280 * qJD(1) - t281 * qJD(2);
	t275 = -qJD(1) * t302 - t294 * t305 + (qJD(2) * t288 + qJD(1)) * t310;
	t274 = -t289 * t300 - t275 * t292 + (-t281 * t289 - t292 * t313) * qJD(5);
	t273 = -t292 * t300 + t275 * t289 + (-t281 * t292 + t289 * t313) * qJD(5);
	t1 = [t295, t276 * t292 - t282 * t304, 0, 0, t273, 0; t274, t278 * t292 - t280 * t304, 0, 0, -t314, 0; 0, (-t290 * t304 + t292 * t305) * t287, 0, 0, -t289 * t299 + (t288 * t289 + t292 * t312) * qJD(5), 0; t314, -t276 * t289 - t282 * t303, 0, 0, -t274, 0; t273, -t278 * t289 - t280 * t303, 0, 0, t295, 0; 0, (-t289 * t305 - t290 * t303) * t287, 0, 0, -t292 * t299 + (t288 * t292 - t289 * t312) * qJD(5), 0; -t278, t275, 0, 0, 0, 0; t276, -t277, 0, 0, 0, 0; 0, -t299, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:43
	% EndTime: 2019-10-10 09:46:44
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (289->73), mult. (898->141), div. (0->0), fcn. (964->10), ass. (0->65)
	t445 = cos(pkin(6));
	t449 = sin(qJ(1));
	t452 = cos(qJ(2));
	t479 = t449 * t452;
	t448 = sin(qJ(2));
	t453 = cos(qJ(1));
	t480 = t448 * t453;
	t435 = t445 * t480 + t479;
	t436 = t445 * t479 + t480;
	t423 = t436 * qJD(1) + t435 * qJD(2);
	t478 = t452 * t453;
	t467 = t445 * t478;
	t482 = t448 * t449;
	t434 = -t467 + t482;
	t447 = sin(qJ(5));
	t451 = cos(qJ(5));
	t444 = sin(pkin(6));
	t483 = t444 * t453;
	t427 = -t434 * t447 + t451 * t483;
	t477 = qJD(1) * t444;
	t466 = t449 * t477;
	t419 = t427 * qJD(5) + t423 * t451 - t447 * t466;
	t462 = qJD(2) * t445 + qJD(1);
	t468 = t445 * t482;
	t476 = qJD(2) * t448;
	t424 = -qJD(1) * t468 - t449 * t476 + t462 * t478;
	t428 = t434 * t451 + t447 * t483;
	t446 = sin(qJ(6));
	t450 = cos(qJ(6));
	t494 = -t419 * t450 - t424 * t446 + (t428 * t446 - t435 * t450) * qJD(6);
	t493 = (t428 * t450 + t435 * t446) * qJD(6) + t419 * t446 - t424 * t450;
	t490 = (qJD(2) * t451 + qJD(6)) * t452;
	t485 = t444 * t449;
	t484 = t444 * t452;
	t481 = t448 * t450;
	t475 = qJD(2) * t452;
	t474 = qJD(5) * t447;
	t473 = qJD(5) * t451;
	t472 = qJD(6) * t446;
	t471 = qJD(6) * t450;
	t470 = qJD(6) * t451;
	t469 = t447 * t485;
	t465 = t453 * t477;
	t464 = t444 * t476;
	t463 = t444 * t475;
	t461 = -qJD(2) - t470;
	t421 = -qJD(1) * t467 - t453 * t475 + t462 * t482;
	t437 = -t468 + t478;
	t459 = -t437 * t470 + t421;
	t458 = -t435 * t470 - t423;
	t430 = -t436 * t447 - t451 * t485;
	t457 = t445 * t447 + t451 * t484;
	t432 = -t445 * t451 + t447 * t484;
	t422 = -t435 * qJD(1) - t436 * qJD(2);
	t455 = qJD(6) * t436 - t422 * t451 + t437 * t474;
	t454 = qJD(6) * t434 - t424 * t451 + t435 * t474;
	t418 = -t428 * qJD(5) - t423 * t447 - t451 * t466;
	t431 = t436 * t451 - t469;
	t426 = t432 * qJD(5) + t451 * t464;
	t425 = t457 * qJD(5) - t447 * t464;
	t417 = t430 * qJD(5) - t421 * t451 - t447 * t465;
	t416 = -t421 * t447 - qJD(5) * t469 + (qJD(5) * t436 + t465) * t451;
	t415 = t417 * t450 + t422 * t446 + (-t431 * t446 + t437 * t450) * qJD(6);
	t414 = -t417 * t446 + t422 * t450 + (-t431 * t450 - t437 * t446) * qJD(6);
	t1 = [t494, t459 * t446 - t455 * t450, 0, 0, -t416 * t450 - t430 * t472, t414; t415, t458 * t446 - t454 * t450, 0, 0, t418 * t450 - t427 * t472, -t493; 0, (t450 * t490 + (t461 * t446 - t450 * t474) * t448) * t444, 0, 0, t425 * t450 - t432 * t472, t450 * t463 - t426 * t446 + (-t444 * t446 * t448 + t450 * t457) * qJD(6); t493, t455 * t446 + t459 * t450, 0, 0, t416 * t446 - t430 * t471, -t415; t414, t454 * t446 + t458 * t450, 0, 0, -t418 * t446 - t427 * t471, t494; 0, (t461 * t481 + (t448 * t474 - t490) * t446) * t444, 0, 0, -t425 * t446 - t432 * t471, -t446 * t463 - t426 * t450 + (-t444 * t481 - t446 * t457) * qJD(6); t418, t422 * t447 + t437 * t473, 0, 0, t417, 0; t416, t424 * t447 + t435 * t473, 0, 0, t419, 0; 0, (t447 * t475 + t448 * t473) * t444, 0, 0, t426, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end