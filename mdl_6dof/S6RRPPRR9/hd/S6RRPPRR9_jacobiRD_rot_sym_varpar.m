% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR9
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
% Datum: 2019-10-10 09:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPRR9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:18
	% EndTime: 2019-10-10 09:50:18
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
	% StartTime: 2019-10-10 09:50:18
	% EndTime: 2019-10-10 09:50:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->12), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t200 = sin(qJ(2));
	t201 = sin(qJ(1));
	t213 = t200 * t201;
	t203 = cos(qJ(1));
	t212 = t200 * t203;
	t202 = cos(qJ(2));
	t211 = t201 * t202;
	t210 = t202 * t203;
	t198 = sin(pkin(6));
	t209 = qJD(1) * t198;
	t208 = qJD(2) * t198;
	t199 = cos(pkin(6));
	t207 = t199 * t210 - t213;
	t206 = -t199 * t211 - t212;
	t205 = -t199 * t212 - t211;
	t204 = t199 * t213 - t210;
	t197 = -t204 * qJD(1) + t207 * qJD(2);
	t196 = t206 * qJD(1) + t205 * qJD(2);
	t195 = t205 * qJD(1) + t206 * qJD(2);
	t194 = -t207 * qJD(1) + t204 * qJD(2);
	t1 = [-t201 * t209, 0, 0, 0, 0, 0; t203 * t209, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t196, t195, 0, 0, 0, 0; -t194, t197, 0, 0, 0, 0; 0, t202 * t208, 0, 0, 0, 0; -t197, t194, 0, 0, 0, 0; t195, t196, 0, 0, 0, 0; 0, -t200 * t208, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:18
	% EndTime: 2019-10-10 09:50:18
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (95->36), mult. (310->71), div. (0->0), fcn. (322->8), ass. (0->37)
	t292 = cos(pkin(6));
	t295 = sin(qJ(1));
	t294 = sin(qJ(2));
	t316 = t295 * t294;
	t307 = t292 * t316;
	t311 = qJD(2) * t294;
	t297 = cos(qJ(2));
	t298 = cos(qJ(1));
	t313 = t298 * t297;
	t284 = -qJD(1) * t307 - t295 * t311 + (qJD(2) * t292 + qJD(1)) * t313;
	t314 = t298 * t294;
	t315 = t295 * t297;
	t286 = t292 * t314 + t315;
	t293 = sin(qJ(5));
	t296 = cos(qJ(5));
	t291 = sin(pkin(6));
	t312 = qJD(1) * t291;
	t306 = t295 * t312;
	t317 = t291 * t298;
	t320 = (-t286 * t293 + t296 * t317) * qJD(5) + t284 * t296 - t293 * t306;
	t319 = t291 * t293;
	t318 = t291 * t296;
	t310 = qJD(5) * t293;
	t309 = qJD(5) * t296;
	t308 = qJD(5) * t297;
	t305 = t298 * t312;
	t304 = t291 * qJD(2) * t297;
	t285 = t292 * t313 - t316;
	t302 = t292 * t315 + t314;
	t301 = t307 - t313;
	t299 = -t296 * t306 - t284 * t293 + (-t286 * t296 - t293 * t317) * qJD(5);
	t283 = t302 * qJD(1) + t286 * qJD(2);
	t282 = t286 * qJD(1) + t302 * qJD(2);
	t281 = -t285 * qJD(1) + t301 * qJD(2);
	t280 = t296 * t305 - t282 * t293 + (-t295 * t319 - t296 * t301) * qJD(5);
	t279 = -t293 * t305 - t282 * t296 + (t293 * t301 - t295 * t318) * qJD(5);
	t1 = [t299, t281 * t293 - t302 * t309, 0, 0, t279, 0; t280, -t283 * t293 + t285 * t309, 0, 0, t320, 0; 0, (-t293 * t311 + t296 * t308) * t291, 0, 0, t296 * t304 + (-t292 * t296 - t294 * t319) * qJD(5), 0; -t320, t281 * t296 + t302 * t310, 0, 0, -t280, 0; t279, -t283 * t296 - t285 * t310, 0, 0, t299, 0; 0, (-t293 * t308 - t296 * t311) * t291, 0, 0, -t293 * t304 + (t292 * t293 - t294 * t318) * qJD(5), 0; t283, t282, 0, 0, 0, 0; t281, -t284, 0, 0, 0, 0; 0, -t304, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:21
	% EndTime: 2019-10-10 09:50:21
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (289->71), mult. (898->138), div. (0->0), fcn. (964->10), ass. (0->63)
	t474 = cos(pkin(6));
	t478 = sin(qJ(1));
	t477 = sin(qJ(2));
	t510 = t477 * t478;
	t497 = t474 * t510;
	t505 = qJD(2) * t477;
	t481 = cos(qJ(2));
	t482 = cos(qJ(1));
	t507 = t481 * t482;
	t454 = -qJD(1) * t497 - t478 * t505 + (qJD(2) * t474 + qJD(1)) * t507;
	t508 = t478 * t481;
	t509 = t477 * t482;
	t465 = t474 * t509 + t508;
	t476 = sin(qJ(5));
	t480 = cos(qJ(5));
	t473 = sin(pkin(6));
	t511 = t473 * t482;
	t459 = t465 * t480 + t476 * t511;
	t506 = qJD(1) * t473;
	t496 = t478 * t506;
	t446 = t459 * qJD(5) + t454 * t476 + t480 * t496;
	t487 = t474 * t508 + t509;
	t453 = t487 * qJD(1) + t465 * qJD(2);
	t498 = t480 * t511;
	t460 = -t465 * t476 + t498;
	t464 = -t474 * t507 + t510;
	t475 = sin(qJ(6));
	t479 = cos(qJ(6));
	t523 = -t446 * t479 + t453 * t475 + (-t460 * t475 + t464 * t479) * qJD(6);
	t522 = (t460 * t479 + t464 * t475) * qJD(6) - t446 * t475 - t453 * t479;
	t502 = qJD(5) * t481;
	t519 = (qJD(2) * t476 + qJD(6)) * t477 - t480 * t502;
	t514 = t473 * t476;
	t513 = t473 * t480;
	t512 = t473 * t481;
	t504 = qJD(5) * t476;
	t503 = qJD(5) * t480;
	t501 = qJD(6) * t475;
	t500 = qJD(6) * t476;
	t499 = qJD(6) * t479;
	t495 = t482 * t506;
	t494 = t473 * t505;
	t493 = qJD(2) * t512;
	t452 = t465 * qJD(1) + t487 * qJD(2);
	t490 = t487 * t500 + t452;
	t489 = t464 * t500 - t454;
	t488 = (-qJD(2) - t500) * t481;
	t486 = t497 - t507;
	t458 = -t476 * t486 + t478 * t513;
	t457 = -t478 * t514 - t480 * t486;
	t462 = -t474 * t476 + t477 * t513;
	t463 = t474 * t480 + t477 * t514;
	t451 = t464 * qJD(1) + t486 * qJD(2);
	t484 = qJD(6) * t486 + t451 * t476 - t487 * t503;
	t483 = qJD(6) * t465 + t453 * t476 + t464 * t503;
	t445 = t454 * t480 + qJD(5) * t498 + (-qJD(5) * t465 - t496) * t476;
	t456 = -t463 * qJD(5) + t480 * t493;
	t455 = t462 * qJD(5) + t476 * t493;
	t449 = t457 * qJD(5) - t452 * t476 + t480 * t495;
	t448 = t458 * qJD(5) + t452 * t480 + t476 * t495;
	t444 = t449 * t479 + t451 * t475 + (-t458 * t475 - t479 * t487) * qJD(6);
	t443 = -t449 * t475 + t451 * t479 + (-t458 * t479 + t475 * t487) * qJD(6);
	t1 = [t523, t490 * t475 + t484 * t479, 0, 0, -t448 * t479 - t457 * t501, t443; t444, t489 * t475 - t483 * t479, 0, 0, t445 * t479 - t459 * t501, t522; 0, (t475 * t488 - t519 * t479) * t473, 0, 0, t456 * t479 - t462 * t501, -t479 * t494 - t455 * t475 + (-t463 * t479 - t475 * t512) * qJD(6); -t522, -t484 * t475 + t490 * t479, 0, 0, t448 * t475 - t457 * t499, -t444; t443, t483 * t475 + t489 * t479, 0, 0, -t445 * t475 - t459 * t499, t523; 0, (t519 * t475 + t479 * t488) * t473, 0, 0, -t456 * t475 - t462 * t499, t475 * t494 - t455 * t479 + (t463 * t475 - t479 * t512) * qJD(6); t445, -t451 * t480 - t487 * t504, 0, 0, t449, 0; t448, t453 * t480 - t464 * t504, 0, 0, t446, 0; 0, (t476 * t502 + t480 * t505) * t473, 0, 0, t455, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end