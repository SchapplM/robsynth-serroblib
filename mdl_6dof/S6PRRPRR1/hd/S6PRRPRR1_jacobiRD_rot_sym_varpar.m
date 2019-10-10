% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(11));
	t58 = sin(pkin(11));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:36
	% EndTime: 2019-10-09 22:25:36
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t211 = sin(pkin(6));
	t214 = sin(qJ(3));
	t227 = t211 * t214;
	t216 = cos(qJ(3));
	t226 = t211 * t216;
	t213 = cos(pkin(6));
	t215 = sin(qJ(2));
	t225 = t213 * t215;
	t217 = cos(qJ(2));
	t224 = t213 * t217;
	t223 = qJD(2) * t215;
	t222 = qJD(3) * t214;
	t221 = qJD(3) * t216;
	t220 = qJD(3) * t217;
	t219 = t211 * qJD(2) * t217;
	t210 = sin(pkin(11));
	t212 = cos(pkin(11));
	t206 = -t210 * t215 + t212 * t224;
	t207 = t210 * t217 + t212 * t225;
	t208 = -t210 * t224 - t212 * t215;
	t218 = t210 * t225 - t212 * t217;
	t205 = t218 * qJD(2);
	t204 = t208 * qJD(2);
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t1 = [0, t205 * t216 - t208 * t222, -t204 * t214 + (-t210 * t227 + t216 * t218) * qJD(3), 0, 0, 0; 0, -t203 * t216 - t206 * t222, -t202 * t214 + (-t207 * t216 + t212 * t227) * qJD(3), 0, 0, 0; 0, (-t214 * t220 - t216 * t223) * t211, -t214 * t219 + (-t213 * t214 - t215 * t226) * qJD(3), 0, 0, 0; 0, -t205 * t214 - t208 * t221, -t204 * t216 + (-t210 * t226 - t214 * t218) * qJD(3), 0, 0, 0; 0, t203 * t214 - t206 * t221, -t202 * t216 + (t207 * t214 + t212 * t226) * qJD(3), 0, 0, 0; 0, (t214 * t223 - t216 * t220) * t211, -t216 * t219 + (-t213 * t216 + t215 * t227) * qJD(3), 0, 0, 0; 0, t204, 0, 0, 0, 0; 0, t202, 0, 0, 0, 0; 0, t219, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:36
	% EndTime: 2019-10-09 22:25:36
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (66->23), mult. (140->60), div. (0->0), fcn. (148->8), ass. (0->28)
	t230 = sin(pkin(11));
	t231 = sin(pkin(6));
	t246 = t230 * t231;
	t232 = cos(pkin(11));
	t245 = t231 * t232;
	t234 = sin(qJ(2));
	t244 = t231 * t234;
	t233 = cos(pkin(6));
	t243 = t233 * t234;
	t235 = cos(qJ(2));
	t242 = t233 * t235;
	t241 = qJD(2) * t234;
	t229 = qJ(3) + pkin(12);
	t227 = sin(t229);
	t240 = qJD(3) * t227;
	t228 = cos(t229);
	t239 = qJD(3) * t228;
	t238 = qJD(3) * t235;
	t237 = t231 * qJD(2) * t235;
	t223 = -t230 * t234 + t232 * t242;
	t224 = t230 * t235 + t232 * t243;
	t225 = -t230 * t242 - t232 * t234;
	t236 = t230 * t243 - t232 * t235;
	t222 = t236 * qJD(2);
	t221 = t225 * qJD(2);
	t220 = t224 * qJD(2);
	t219 = t223 * qJD(2);
	t1 = [0, t222 * t228 - t225 * t240, -t221 * t227 + (-t227 * t246 + t228 * t236) * qJD(3), 0, 0, 0; 0, -t220 * t228 - t223 * t240, -t219 * t227 + (-t224 * t228 + t227 * t245) * qJD(3), 0, 0, 0; 0, (-t227 * t238 - t228 * t241) * t231, -t227 * t237 + (-t227 * t233 - t228 * t244) * qJD(3), 0, 0, 0; 0, -t222 * t227 - t225 * t239, -t221 * t228 + (-t227 * t236 - t228 * t246) * qJD(3), 0, 0, 0; 0, t220 * t227 - t223 * t239, -t219 * t228 + (t224 * t227 + t228 * t245) * qJD(3), 0, 0, 0; 0, (t227 * t241 - t228 * t238) * t231, -t228 * t237 + (t227 * t244 - t228 * t233) * qJD(3), 0, 0, 0; 0, t221, 0, 0, 0, 0; 0, t219, 0, 0, 0, 0; 0, t237, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:36
	% EndTime: 2019-10-09 22:25:36
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (182->21), mult. (212->50), div. (0->0), fcn. (224->8), ass. (0->37)
	t283 = qJ(3) + pkin(12) + qJ(5);
	t281 = sin(t283);
	t284 = qJD(3) + qJD(5);
	t303 = t281 * t284;
	t282 = cos(t283);
	t302 = t282 * t284;
	t286 = sin(pkin(6));
	t301 = t284 * t286;
	t290 = cos(qJ(2));
	t300 = t284 * t290;
	t288 = cos(pkin(6));
	t289 = sin(qJ(2));
	t299 = t288 * t289;
	t298 = t288 * t290;
	t297 = qJD(2) * t289;
	t296 = t289 * t301;
	t295 = t286 * qJD(2) * t290;
	t285 = sin(pkin(11));
	t287 = cos(pkin(11));
	t277 = -t285 * t289 + t287 * t298;
	t273 = t277 * qJD(2);
	t294 = t287 * t301 - t273;
	t279 = -t285 * t298 - t287 * t289;
	t275 = t279 * qJD(2);
	t293 = -t285 * t301 - t275;
	t278 = t285 * t290 + t287 * t299;
	t292 = t285 * t299 - t287 * t290;
	t291 = -t284 * t288 - t295;
	t276 = t292 * qJD(2);
	t274 = t278 * qJD(2);
	t272 = t281 * t296 + t291 * t282;
	t271 = t291 * t281 - t282 * t296;
	t270 = t293 * t282 - t292 * t303;
	t269 = t293 * t281 + t292 * t302;
	t268 = t278 * t303 + t294 * t282;
	t267 = -t278 * t302 + t294 * t281;
	t1 = [0, t276 * t282 - t279 * t303, t269, 0, t269, 0; 0, -t274 * t282 - t277 * t303, t267, 0, t267, 0; 0, (-t281 * t300 - t282 * t297) * t286, t271, 0, t271, 0; 0, -t276 * t281 - t279 * t302, t270, 0, t270, 0; 0, t274 * t281 - t277 * t302, t268, 0, t268, 0; 0, (t281 * t297 - t282 * t300) * t286, t272, 0, t272, 0; 0, t275, 0, 0, 0, 0; 0, t273, 0, 0, 0, 0; 0, t295, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:38
	% EndTime: 2019-10-09 22:25:39
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (520->64), mult. (672->131), div. (0->0), fcn. (736->10), ass. (0->66)
	t453 = qJ(3) + pkin(12) + qJ(5);
	t451 = sin(t453);
	t452 = cos(t453);
	t460 = sin(qJ(2));
	t454 = qJD(3) + qJD(5);
	t462 = cos(qJ(2));
	t489 = t454 * t462;
	t492 = (qJD(2) * t452 - qJD(6)) * t460 + t451 * t489;
	t491 = t451 * t454;
	t490 = t452 * t454;
	t455 = sin(pkin(11));
	t456 = sin(pkin(6));
	t488 = t455 * t456;
	t457 = cos(pkin(11));
	t487 = t456 * t457;
	t486 = t456 * t460;
	t485 = t456 * t462;
	t458 = cos(pkin(6));
	t484 = t458 * t460;
	t483 = t458 * t462;
	t482 = qJD(2) * t460;
	t481 = qJD(2) * t462;
	t480 = qJD(6) * t452;
	t459 = sin(qJ(6));
	t479 = qJD(6) * t459;
	t461 = cos(qJ(6));
	t478 = qJD(6) * t461;
	t476 = t455 * t484;
	t475 = t451 * t486;
	t474 = t452 * t486;
	t473 = t456 * t482;
	t466 = -t455 * t460 + t457 * t483;
	t441 = t466 * qJD(2);
	t471 = t454 * t487 - t441;
	t447 = t455 * t483 + t457 * t460;
	t443 = t447 * qJD(2);
	t470 = t454 * t488 - t443;
	t469 = -t466 * t480 + t441;
	t468 = t447 * t480 - t443;
	t467 = (qJD(2) - t480) * t462;
	t446 = t455 * t462 + t457 * t484;
	t465 = t454 * t458 + t456 * t481;
	t442 = t446 * qJD(2);
	t464 = qJD(6) * t446 - t442 * t452 - t466 * t491;
	t444 = -qJD(2) * t476 + t457 * t481;
	t448 = t457 * t462 - t476;
	t463 = qJD(6) * t448 - t444 * t452 + t447 * t491;
	t440 = t458 * t451 + t474;
	t439 = t458 * t452 - t475;
	t438 = t448 * t452 + t451 * t488;
	t437 = -t448 * t451 + t452 * t488;
	t436 = t446 * t452 - t451 * t487;
	t435 = -t446 * t451 - t452 * t487;
	t434 = t465 * t452 - t454 * t475;
	t433 = -t465 * t451 - t454 * t474;
	t432 = -t448 * t491 + t470 * t452;
	t431 = -t448 * t490 - t470 * t451;
	t430 = -t446 * t491 - t471 * t452;
	t429 = -t446 * t490 + t471 * t451;
	t428 = t433 * t461 - t439 * t479;
	t427 = -t433 * t459 - t439 * t478;
	t426 = t431 * t461 - t437 * t479;
	t425 = -t431 * t459 - t437 * t478;
	t424 = t429 * t461 - t435 * t479;
	t423 = -t429 * t459 - t435 * t478;
	t1 = [0, t468 * t459 + t463 * t461, t426, 0, t426, -t432 * t459 + t444 * t461 + (-t438 * t461 - t447 * t459) * qJD(6); 0, t469 * t459 + t464 * t461, t424, 0, t424, -t430 * t459 + t442 * t461 + (-t436 * t461 + t459 * t466) * qJD(6); 0, (t459 * t467 - t492 * t461) * t456, t428, 0, t428, t461 * t473 - t434 * t459 + (-t440 * t461 + t459 * t485) * qJD(6); 0, -t463 * t459 + t468 * t461, t425, 0, t425, -t432 * t461 - t444 * t459 + (t438 * t459 - t447 * t461) * qJD(6); 0, -t464 * t459 + t469 * t461, t423, 0, t423, -t430 * t461 - t442 * t459 + (t436 * t459 + t461 * t466) * qJD(6); 0, (t492 * t459 + t461 * t467) * t456, t427, 0, t427, -t459 * t473 - t434 * t461 + (t440 * t459 + t461 * t485) * qJD(6); 0, -t444 * t451 - t447 * t490, t432, 0, t432, 0; 0, -t442 * t451 + t466 * t490, t430, 0, t430, 0; 0, (-t451 * t482 + t452 * t489) * t456, t434, 0, t434, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end