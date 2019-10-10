% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR2
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
% Datum: 2019-10-09 22:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
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
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
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
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:35
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
	% StartTime: 2019-10-09 22:48:36
	% EndTime: 2019-10-09 22:48:36
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (214->34), mult. (354->74), div. (0->0), fcn. (372->10), ass. (0->47)
	t388 = sin(pkin(11));
	t391 = cos(pkin(11));
	t394 = cos(qJ(2));
	t392 = cos(pkin(6));
	t393 = sin(qJ(2));
	t404 = t392 * t393;
	t380 = t388 * t394 + t391 * t404;
	t386 = qJ(3) + qJ(4);
	t383 = sin(t386);
	t403 = t392 * t394;
	t379 = -t388 * t393 + t391 * t403;
	t375 = t379 * qJD(2);
	t385 = qJD(3) + qJD(4);
	t389 = sin(pkin(6));
	t406 = t385 * t389;
	t400 = t391 * t406 - t375;
	t384 = cos(t386);
	t408 = t384 * t385;
	t368 = -t380 * t408 + t400 * t383;
	t387 = sin(pkin(12));
	t412 = t368 * t387;
	t396 = t388 * t404 - t391 * t394;
	t381 = -t388 * t403 - t391 * t393;
	t377 = t381 * qJD(2);
	t399 = t388 * t406 + t377;
	t370 = -t399 * t383 + t396 * t408;
	t411 = t370 * t387;
	t395 = qJD(2) * t389 * t394 + t385 * t392;
	t401 = t393 * t406;
	t373 = -t395 * t383 - t384 * t401;
	t410 = t373 * t387;
	t409 = t383 * t385;
	t407 = t384 * t393;
	t405 = t385 * t394;
	t402 = t383 * t405;
	t376 = t380 * qJD(2);
	t398 = t376 * t384 + t379 * t409;
	t378 = t396 * qJD(2);
	t397 = -t378 * t384 + t381 * t409;
	t390 = cos(pkin(12));
	t374 = -t383 * t401 + t395 * t384;
	t372 = t373 * t390;
	t371 = t399 * t384 + t396 * t409;
	t369 = -t380 * t409 - t400 * t384;
	t367 = t370 * t390;
	t366 = t368 * t390;
	t1 = [0, t377 * t387 - t397 * t390, t367, t367, 0, 0; 0, t375 * t387 - t398 * t390, t366, t366, 0, 0; 0, (-t390 * t402 + (t387 * t394 - t390 * t407) * qJD(2)) * t389, t372, t372, 0, 0; 0, t377 * t390 + t397 * t387, -t411, -t411, 0, 0; 0, t375 * t390 + t398 * t387, -t412, -t412, 0, 0; 0, (t387 * t402 + (t387 * t407 + t390 * t394) * qJD(2)) * t389, -t410, -t410, 0, 0; 0, t378 * t383 + t381 * t408, t371, t371, 0, 0; 0, -t376 * t383 + t379 * t408, t369, t369, 0, 0; 0, (-qJD(2) * t383 * t393 + t384 * t405) * t389, t374, t374, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:37
	% EndTime: 2019-10-09 22:48:37
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (466->65), mult. (672->131), div. (0->0), fcn. (736->10), ass. (0->67)
	t459 = qJ(3) + qJ(4);
	t455 = sin(t459);
	t456 = cos(t459);
	t464 = sin(qJ(2));
	t458 = qJD(3) + qJD(4);
	t465 = cos(qJ(2));
	t492 = t458 * t465;
	t495 = (qJD(2) * t456 - qJD(6)) * t464 + t455 * t492;
	t494 = t455 * t458;
	t493 = t456 * t458;
	t460 = sin(pkin(11));
	t461 = sin(pkin(6));
	t491 = t460 * t461;
	t462 = cos(pkin(11));
	t490 = t461 * t462;
	t489 = t461 * t464;
	t488 = t461 * t465;
	t463 = cos(pkin(6));
	t487 = t463 * t464;
	t486 = t463 * t465;
	t485 = qJD(2) * t464;
	t484 = qJD(2) * t465;
	t457 = pkin(12) + qJ(6);
	t453 = sin(t457);
	t483 = qJD(6) * t453;
	t454 = cos(t457);
	t482 = qJD(6) * t454;
	t481 = qJD(6) * t456;
	t479 = t460 * t487;
	t478 = t455 * t489;
	t477 = t456 * t489;
	t476 = t461 * t485;
	t469 = -t460 * t464 + t462 * t486;
	t443 = t469 * qJD(2);
	t474 = t458 * t490 - t443;
	t449 = t460 * t486 + t462 * t464;
	t445 = t449 * qJD(2);
	t473 = t458 * t491 - t445;
	t472 = -t469 * t481 + t443;
	t471 = t449 * t481 - t445;
	t470 = (qJD(2) - t481) * t465;
	t448 = t460 * t465 + t462 * t487;
	t468 = t458 * t463 + t461 * t484;
	t444 = t448 * qJD(2);
	t467 = qJD(6) * t448 - t444 * t456 - t469 * t494;
	t446 = -qJD(2) * t479 + t462 * t484;
	t450 = t462 * t465 - t479;
	t466 = qJD(6) * t450 - t446 * t456 + t449 * t494;
	t442 = t463 * t455 + t477;
	t441 = t463 * t456 - t478;
	t440 = t450 * t456 + t455 * t491;
	t439 = -t450 * t455 + t456 * t491;
	t438 = t448 * t456 - t455 * t490;
	t437 = -t448 * t455 - t456 * t490;
	t436 = t468 * t456 - t458 * t478;
	t435 = -t468 * t455 - t458 * t477;
	t434 = -t450 * t494 + t473 * t456;
	t433 = -t450 * t493 - t473 * t455;
	t432 = -t448 * t494 - t474 * t456;
	t431 = -t448 * t493 + t474 * t455;
	t430 = t435 * t454 - t441 * t483;
	t429 = -t435 * t453 - t441 * t482;
	t428 = t433 * t454 - t439 * t483;
	t427 = -t433 * t453 - t439 * t482;
	t426 = t431 * t454 - t437 * t483;
	t425 = -t431 * t453 - t437 * t482;
	t1 = [0, t471 * t453 + t466 * t454, t428, t428, 0, -t434 * t453 + t446 * t454 + (-t440 * t454 - t449 * t453) * qJD(6); 0, t472 * t453 + t467 * t454, t426, t426, 0, -t432 * t453 + t444 * t454 + (-t438 * t454 + t453 * t469) * qJD(6); 0, (t453 * t470 - t495 * t454) * t461, t430, t430, 0, t454 * t476 - t436 * t453 + (-t442 * t454 + t453 * t488) * qJD(6); 0, -t466 * t453 + t471 * t454, t427, t427, 0, -t434 * t454 - t446 * t453 + (t440 * t453 - t449 * t454) * qJD(6); 0, -t467 * t453 + t472 * t454, t425, t425, 0, -t432 * t454 - t444 * t453 + (t438 * t453 + t454 * t469) * qJD(6); 0, (t495 * t453 + t454 * t470) * t461, t429, t429, 0, -t453 * t476 - t436 * t454 + (t442 * t453 + t454 * t488) * qJD(6); 0, -t446 * t455 - t449 * t493, t434, t434, 0, 0; 0, -t444 * t455 + t469 * t493, t432, t432, 0, 0; 0, (-t455 * t485 + t456 * t492) * t461, t436, t436, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end