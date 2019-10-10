% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
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
	% StartTime: 2019-10-10 10:24:32
	% EndTime: 2019-10-10 10:24:32
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
	% StartTime: 2019-10-10 10:24:32
	% EndTime: 2019-10-10 10:24:32
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
	t280 = t284 * qJD(1) + t283 * qJD(2);
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
	t1 = [t297, t279 * t291 + t285 * t305, 0, t276, 0, 0; t277, t281 * t291 + t283 * t305, 0, t316, 0, 0; 0, (t291 * t307 + t292 * t305) * t289, 0, t294 * t301 + (-t290 * t294 + t291 * t314) * qJD(4), 0, 0; -t316, t279 * t294 - t285 * t306, 0, -t277, 0, 0; t276, t281 * t294 - t283 * t306, 0, t297, 0, 0; 0, (-t292 * t306 + t294 * t307) * t289, 0, -t291 * t301 + (t290 * t291 + t294 * t314) * qJD(4), 0, 0; -t281, t278, 0, 0, 0, 0; t279, -t280, 0, 0, 0, 0; 0, -t301, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:32
	% EndTime: 2019-10-10 10:24:32
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (145->38), mult. (310->71), div. (0->0), fcn. (322->8), ass. (0->38)
	t310 = cos(pkin(6));
	t312 = sin(qJ(1));
	t313 = cos(qJ(2));
	t329 = t312 * t313;
	t311 = sin(qJ(2));
	t314 = cos(qJ(1));
	t330 = t311 * t314;
	t300 = t310 * t330 + t329;
	t301 = t310 * t329 + t330;
	t297 = t301 * qJD(1) + t300 * qJD(2);
	t308 = qJ(4) + pkin(11);
	t306 = sin(t308);
	t307 = cos(t308);
	t328 = t313 * t314;
	t322 = t310 * t328;
	t331 = t311 * t312;
	t317 = t322 - t331;
	t309 = sin(pkin(6));
	t327 = qJD(1) * t309;
	t321 = t312 * t327;
	t332 = t309 * t314;
	t335 = (t306 * t317 + t307 * t332) * qJD(4) + t297 * t307 - t306 * t321;
	t334 = t309 * t312;
	t333 = t309 * t313;
	t326 = qJD(2) * t313;
	t325 = qJD(4) * t306;
	t324 = qJD(4) * t307;
	t323 = qJD(4) * t311;
	t320 = t314 * t327;
	t319 = t309 * qJD(2) * t311;
	t302 = -t310 * t331 + t328;
	t315 = -t307 * t321 - t297 * t306 + (-t306 * t332 + t307 * t317) * qJD(4);
	t298 = t302 * qJD(1) + t317 * qJD(2);
	t296 = -t300 * qJD(1) - t301 * qJD(2);
	t295 = -qJD(1) * t322 - t314 * t326 + (qJD(2) * t310 + qJD(1)) * t331;
	t294 = t307 * t320 - t295 * t306 + (t301 * t307 - t306 * t334) * qJD(4);
	t293 = -t306 * t320 - t295 * t307 + (-t301 * t306 - t307 * t334) * qJD(4);
	t1 = [t315, t296 * t306 + t302 * t324, 0, t293, 0, 0; t294, t298 * t306 + t300 * t324, 0, t335, 0, 0; 0, (t306 * t326 + t307 * t323) * t309, 0, t307 * t319 + (t306 * t333 - t307 * t310) * qJD(4), 0, 0; -t335, t296 * t307 - t302 * t325, 0, -t294, 0, 0; t293, t298 * t307 - t300 * t325, 0, t315, 0, 0; 0, (-t306 * t323 + t307 * t326) * t309, 0, -t306 * t319 + (t306 * t310 + t307 * t333) * qJD(4), 0, 0; -t298, t295, 0, 0, 0, 0; t296, -t297, 0, 0, 0, 0; 0, -t319, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:35
	% EndTime: 2019-10-10 10:24:35
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (424->73), mult. (898->141), div. (0->0), fcn. (964->10), ass. (0->65)
	t485 = cos(pkin(6));
	t487 = sin(qJ(2));
	t491 = cos(qJ(1));
	t517 = t491 * t487;
	t488 = sin(qJ(1));
	t490 = cos(qJ(2));
	t518 = t488 * t490;
	t474 = t485 * t517 + t518;
	t475 = t485 * t518 + t517;
	t462 = t475 * qJD(1) + t474 * qJD(2);
	t516 = t491 * t490;
	t519 = t488 * t487;
	t473 = -t485 * t516 + t519;
	t483 = qJ(4) + pkin(11);
	t481 = sin(t483);
	t482 = cos(t483);
	t484 = sin(pkin(6));
	t521 = t484 * t491;
	t468 = t473 * t482 + t481 * t521;
	t515 = qJD(1) * t484;
	t504 = t488 * t515;
	t455 = t468 * qJD(4) + t462 * t481 + t482 * t504;
	t506 = t485 * t519;
	t514 = qJD(2) * t487;
	t463 = -qJD(1) * t506 - t488 * t514 + (qJD(2) * t485 + qJD(1)) * t516;
	t505 = t482 * t521;
	t469 = -t473 * t481 + t505;
	t486 = sin(qJ(6));
	t489 = cos(qJ(6));
	t532 = -t455 * t489 + (-t469 * t486 - t474 * t489) * qJD(6) - t463 * t486;
	t531 = (t469 * t489 - t474 * t486) * qJD(6) - t455 * t486 + t463 * t489;
	t528 = (qJD(2) * t481 + qJD(6)) * t490;
	t523 = t484 * t488;
	t522 = t484 * t490;
	t520 = t487 * t489;
	t513 = qJD(2) * t490;
	t512 = qJD(4) * t481;
	t511 = qJD(4) * t482;
	t510 = qJD(4) * t487;
	t509 = qJD(6) * t481;
	t508 = qJD(6) * t486;
	t507 = qJD(6) * t489;
	t503 = t491 * t515;
	t502 = t484 * t514;
	t501 = t484 * t513;
	t500 = -qJD(2) - t509;
	t495 = t506 - t516;
	t460 = t473 * qJD(1) + t495 * qJD(2);
	t498 = t495 * t509 + t460;
	t497 = -t474 * t509 - t462;
	t467 = t475 * t481 + t482 * t523;
	t466 = t475 * t482 - t481 * t523;
	t471 = -t485 * t481 - t482 * t522;
	t496 = t481 * t522 - t485 * t482;
	t461 = -t474 * qJD(1) - t475 * qJD(2);
	t493 = -qJD(6) * t475 + t461 * t481 - t495 * t511;
	t492 = -qJD(6) * t473 + t463 * t481 + t474 * t511;
	t454 = t462 * t482 + qJD(4) * t505 + (-qJD(4) * t473 - t504) * t481;
	t465 = t496 * qJD(4) + t482 * t502;
	t464 = t471 * qJD(4) + t481 * t502;
	t458 = t466 * qJD(4) - t460 * t481 + t482 * t503;
	t457 = t467 * qJD(4) + t460 * t482 + t481 * t503;
	t453 = t458 * t489 + t461 * t486 + (-t467 * t486 - t489 * t495) * qJD(6);
	t452 = -t458 * t486 + t461 * t489 + (-t467 * t489 + t486 * t495) * qJD(6);
	t1 = [t532, t498 * t486 + t493 * t489, 0, -t457 * t489 - t466 * t508, 0, t452; t453, t497 * t486 + t492 * t489, 0, t454 * t489 - t468 * t508, 0, t531; 0, (t489 * t528 + (t500 * t486 + t489 * t511) * t487) * t484, 0, t465 * t489 - t471 * t508, 0, t489 * t501 - t464 * t486 + (-t484 * t486 * t487 + t489 * t496) * qJD(6); -t531, -t493 * t486 + t498 * t489, 0, t457 * t486 - t466 * t507, 0, -t453; t452, -t492 * t486 + t497 * t489, 0, -t454 * t486 - t468 * t507, 0, t532; 0, (t500 * t520 + (-t482 * t510 - t528) * t486) * t484, 0, -t465 * t486 - t471 * t507, 0, -t486 * t501 - t464 * t489 + (-t484 * t520 - t486 * t496) * qJD(6); t454, -t461 * t482 - t495 * t512, 0, t458, 0, 0; t457, -t463 * t482 + t474 * t512, 0, t455, 0, 0; 0, (t481 * t510 - t482 * t513) * t484, 0, t464, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end