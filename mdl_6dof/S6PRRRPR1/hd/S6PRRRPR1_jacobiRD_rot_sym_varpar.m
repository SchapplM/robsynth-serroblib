% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
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
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:44
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
	% StartTime: 2019-10-09 22:46:44
	% EndTime: 2019-10-09 22:46:44
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (134->21), mult. (212->50), div. (0->0), fcn. (224->8), ass. (0->37)
	t274 = qJ(3) + qJ(4);
	t271 = sin(t274);
	t273 = qJD(3) + qJD(4);
	t293 = t271 * t273;
	t272 = cos(t274);
	t292 = t272 * t273;
	t276 = sin(pkin(6));
	t291 = t273 * t276;
	t280 = cos(qJ(2));
	t290 = t273 * t280;
	t278 = cos(pkin(6));
	t279 = sin(qJ(2));
	t289 = t278 * t279;
	t288 = t278 * t280;
	t287 = qJD(2) * t279;
	t286 = t279 * t291;
	t285 = t276 * qJD(2) * t280;
	t275 = sin(pkin(11));
	t277 = cos(pkin(11));
	t267 = -t275 * t279 + t277 * t288;
	t263 = t267 * qJD(2);
	t284 = t277 * t291 - t263;
	t269 = -t275 * t288 - t277 * t279;
	t265 = t269 * qJD(2);
	t283 = -t275 * t291 - t265;
	t268 = t275 * t280 + t277 * t289;
	t282 = t275 * t289 - t277 * t280;
	t281 = -t273 * t278 - t285;
	t266 = t282 * qJD(2);
	t264 = t268 * qJD(2);
	t262 = t271 * t286 + t281 * t272;
	t261 = t281 * t271 - t272 * t286;
	t260 = t283 * t272 - t282 * t293;
	t259 = t283 * t271 + t282 * t292;
	t258 = t268 * t293 + t284 * t272;
	t257 = -t268 * t292 + t284 * t271;
	t1 = [0, t266 * t272 - t269 * t293, t259, t259, 0, 0; 0, -t264 * t272 - t267 * t293, t257, t257, 0, 0; 0, (-t271 * t290 - t272 * t287) * t276, t261, t261, 0, 0; 0, -t266 * t271 - t269 * t292, t260, t260, 0, 0; 0, t264 * t271 - t267 * t292, t258, t258, 0, 0; 0, (t271 * t287 - t272 * t290) * t276, t262, t262, 0, 0; 0, t265, 0, 0, 0, 0; 0, t263, 0, 0, 0, 0; 0, t285, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:44
	% EndTime: 2019-10-09 22:46:44
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (182->21), mult. (212->50), div. (0->0), fcn. (224->8), ass. (0->37)
	t285 = qJ(3) + qJ(4) + pkin(12);
	t283 = sin(t285);
	t286 = qJD(3) + qJD(4);
	t305 = t283 * t286;
	t284 = cos(t285);
	t304 = t284 * t286;
	t288 = sin(pkin(6));
	t303 = t286 * t288;
	t292 = cos(qJ(2));
	t302 = t286 * t292;
	t290 = cos(pkin(6));
	t291 = sin(qJ(2));
	t301 = t290 * t291;
	t300 = t290 * t292;
	t299 = qJD(2) * t291;
	t298 = t291 * t303;
	t297 = t288 * qJD(2) * t292;
	t287 = sin(pkin(11));
	t289 = cos(pkin(11));
	t279 = -t287 * t291 + t289 * t300;
	t275 = t279 * qJD(2);
	t296 = t289 * t303 - t275;
	t281 = -t287 * t300 - t289 * t291;
	t277 = t281 * qJD(2);
	t295 = -t287 * t303 - t277;
	t280 = t287 * t292 + t289 * t301;
	t294 = t287 * t301 - t289 * t292;
	t293 = -t286 * t290 - t297;
	t278 = t294 * qJD(2);
	t276 = t280 * qJD(2);
	t274 = t283 * t298 + t293 * t284;
	t273 = t293 * t283 - t284 * t298;
	t272 = t295 * t284 - t294 * t305;
	t271 = t295 * t283 + t294 * t304;
	t270 = t280 * t305 + t296 * t284;
	t269 = -t280 * t304 + t296 * t283;
	t1 = [0, t278 * t284 - t281 * t305, t271, t271, 0, 0; 0, -t276 * t284 - t279 * t305, t269, t269, 0, 0; 0, (-t283 * t302 - t284 * t299) * t288, t273, t273, 0, 0; 0, -t278 * t283 - t281 * t304, t272, t272, 0, 0; 0, t276 * t283 - t279 * t304, t270, t270, 0, 0; 0, (t283 * t299 - t284 * t302) * t288, t274, t274, 0, 0; 0, t277, 0, 0, 0, 0; 0, t275, 0, 0, 0, 0; 0, t297, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:46
	% EndTime: 2019-10-09 22:46:46
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (520->64), mult. (672->131), div. (0->0), fcn. (736->10), ass. (0->66)
	t456 = qJ(3) + qJ(4) + pkin(12);
	t454 = sin(t456);
	t455 = cos(t456);
	t463 = sin(qJ(2));
	t457 = qJD(3) + qJD(4);
	t465 = cos(qJ(2));
	t492 = t457 * t465;
	t495 = (qJD(2) * t455 - qJD(6)) * t463 + t454 * t492;
	t494 = t454 * t457;
	t493 = t455 * t457;
	t458 = sin(pkin(11));
	t459 = sin(pkin(6));
	t491 = t458 * t459;
	t460 = cos(pkin(11));
	t490 = t459 * t460;
	t489 = t459 * t463;
	t488 = t459 * t465;
	t461 = cos(pkin(6));
	t487 = t461 * t463;
	t486 = t461 * t465;
	t485 = qJD(2) * t463;
	t484 = qJD(2) * t465;
	t483 = qJD(6) * t455;
	t462 = sin(qJ(6));
	t482 = qJD(6) * t462;
	t464 = cos(qJ(6));
	t481 = qJD(6) * t464;
	t479 = t458 * t487;
	t478 = t454 * t489;
	t477 = t455 * t489;
	t476 = t459 * t485;
	t469 = -t458 * t463 + t460 * t486;
	t444 = t469 * qJD(2);
	t474 = t457 * t490 - t444;
	t450 = t458 * t486 + t460 * t463;
	t446 = t450 * qJD(2);
	t473 = t457 * t491 - t446;
	t472 = -t469 * t483 + t444;
	t471 = t450 * t483 - t446;
	t470 = (qJD(2) - t483) * t465;
	t449 = t458 * t465 + t460 * t487;
	t468 = t457 * t461 + t459 * t484;
	t445 = t449 * qJD(2);
	t467 = qJD(6) * t449 - t445 * t455 - t469 * t494;
	t447 = -qJD(2) * t479 + t460 * t484;
	t451 = t460 * t465 - t479;
	t466 = qJD(6) * t451 - t447 * t455 + t450 * t494;
	t443 = t461 * t454 + t477;
	t442 = t461 * t455 - t478;
	t441 = t451 * t455 + t454 * t491;
	t440 = -t451 * t454 + t455 * t491;
	t439 = t449 * t455 - t454 * t490;
	t438 = -t449 * t454 - t455 * t490;
	t437 = t468 * t455 - t457 * t478;
	t436 = -t468 * t454 - t457 * t477;
	t435 = -t451 * t494 + t473 * t455;
	t434 = -t451 * t493 - t473 * t454;
	t433 = -t449 * t494 - t474 * t455;
	t432 = -t449 * t493 + t474 * t454;
	t431 = t436 * t464 - t442 * t482;
	t430 = -t436 * t462 - t442 * t481;
	t429 = t434 * t464 - t440 * t482;
	t428 = -t434 * t462 - t440 * t481;
	t427 = t432 * t464 - t438 * t482;
	t426 = -t432 * t462 - t438 * t481;
	t1 = [0, t471 * t462 + t466 * t464, t429, t429, 0, -t435 * t462 + t447 * t464 + (-t441 * t464 - t450 * t462) * qJD(6); 0, t472 * t462 + t467 * t464, t427, t427, 0, -t433 * t462 + t445 * t464 + (-t439 * t464 + t462 * t469) * qJD(6); 0, (t462 * t470 - t495 * t464) * t459, t431, t431, 0, t464 * t476 - t437 * t462 + (-t443 * t464 + t462 * t488) * qJD(6); 0, -t466 * t462 + t471 * t464, t428, t428, 0, -t435 * t464 - t447 * t462 + (t441 * t462 - t450 * t464) * qJD(6); 0, -t467 * t462 + t472 * t464, t426, t426, 0, -t433 * t464 - t445 * t462 + (t439 * t462 + t464 * t469) * qJD(6); 0, (t495 * t462 + t464 * t470) * t459, t430, t430, 0, -t462 * t476 - t437 * t464 + (t443 * t462 + t464 * t488) * qJD(6); 0, -t447 * t454 - t450 * t493, t435, t435, 0, 0; 0, -t445 * t454 + t469 * t493, t433, t433, 0, 0; 0, (-t454 * t485 + t455 * t492) * t459, t437, t437, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end