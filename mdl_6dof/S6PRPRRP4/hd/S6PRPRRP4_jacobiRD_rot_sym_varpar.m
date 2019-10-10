% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRRP4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(10));
	t58 = sin(pkin(10));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->7), mult. (42->23), div. (0->0), fcn. (42->8), ass. (0->14)
	t176 = cos(pkin(6));
	t177 = sin(qJ(2));
	t182 = t176 * t177;
	t178 = cos(qJ(2));
	t181 = t176 * t178;
	t180 = qJD(2) * sin(pkin(6));
	t179 = t177 * t180;
	t175 = cos(pkin(10));
	t174 = cos(pkin(11));
	t172 = sin(pkin(10));
	t171 = sin(pkin(11));
	t170 = (t172 * t182 - t175 * t178) * qJD(2);
	t169 = (-t172 * t178 - t175 * t182) * qJD(2);
	t1 = [0, t170 * t174, 0, 0, 0, 0; 0, t169 * t174, 0, 0, 0, 0; 0, -t174 * t179, 0, 0, 0, 0; 0, -t170 * t171, 0, 0, 0, 0; 0, -t169 * t171, 0, 0, 0, 0; 0, t171 * t179, 0, 0, 0, 0; 0, (-t172 * t181 - t175 * t177) * qJD(2), 0, 0, 0, 0; 0, (-t172 * t177 + t175 * t181) * qJD(2), 0, 0, 0, 0; 0, t178 * t180, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:19
	% EndTime: 2019-10-09 21:48:19
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (66->23), mult. (140->60), div. (0->0), fcn. (148->8), ass. (0->28)
	t224 = sin(pkin(10));
	t225 = sin(pkin(6));
	t240 = t224 * t225;
	t226 = cos(pkin(10));
	t239 = t225 * t226;
	t228 = sin(qJ(2));
	t238 = t225 * t228;
	t227 = cos(pkin(6));
	t237 = t227 * t228;
	t229 = cos(qJ(2));
	t236 = t227 * t229;
	t235 = qJD(2) * t228;
	t223 = pkin(11) + qJ(4);
	t221 = sin(t223);
	t234 = qJD(4) * t221;
	t222 = cos(t223);
	t233 = qJD(4) * t222;
	t232 = qJD(4) * t229;
	t231 = t225 * qJD(2) * t229;
	t217 = -t224 * t228 + t226 * t236;
	t218 = t224 * t229 + t226 * t237;
	t219 = -t224 * t236 - t226 * t228;
	t230 = t224 * t237 - t226 * t229;
	t216 = t230 * qJD(2);
	t215 = t219 * qJD(2);
	t214 = t218 * qJD(2);
	t213 = t217 * qJD(2);
	t1 = [0, t216 * t222 - t219 * t234, 0, -t215 * t221 + (-t221 * t240 + t222 * t230) * qJD(4), 0, 0; 0, -t214 * t222 - t217 * t234, 0, -t213 * t221 + (-t218 * t222 + t221 * t239) * qJD(4), 0, 0; 0, (-t221 * t232 - t222 * t235) * t225, 0, -t221 * t231 + (-t221 * t227 - t222 * t238) * qJD(4), 0, 0; 0, -t216 * t221 - t219 * t233, 0, -t215 * t222 + (-t221 * t230 - t222 * t240) * qJD(4), 0, 0; 0, t214 * t221 - t217 * t233, 0, -t213 * t222 + (t218 * t221 + t222 * t239) * qJD(4), 0, 0; 0, (t221 * t235 - t222 * t232) * t225, 0, -t222 * t231 + (t221 * t238 - t222 * t227) * qJD(4), 0, 0; 0, t215, 0, 0, 0, 0; 0, t213, 0, 0, 0, 0; 0, t231, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:20
	% EndTime: 2019-10-09 21:48:20
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (246->60), mult. (516->128), div. (0->0), fcn. (564->10), ass. (0->55)
	t367 = pkin(11) + qJ(4);
	t365 = sin(t367);
	t366 = cos(t367);
	t373 = sin(qJ(2));
	t375 = cos(qJ(2));
	t391 = qJD(4) * t375;
	t402 = (qJD(2) * t366 - qJD(5)) * t373 + t365 * t391;
	t368 = sin(pkin(10));
	t369 = sin(pkin(6));
	t401 = t368 * t369;
	t370 = cos(pkin(10));
	t400 = t369 * t370;
	t399 = t369 * t373;
	t398 = t369 * t375;
	t371 = cos(pkin(6));
	t397 = t371 * t373;
	t396 = t371 * t375;
	t395 = qJD(2) * t373;
	t394 = qJD(2) * t375;
	t393 = qJD(4) * t365;
	t392 = qJD(4) * t366;
	t390 = qJD(5) * t366;
	t372 = sin(qJ(5));
	t389 = qJD(5) * t372;
	t374 = cos(qJ(5));
	t388 = qJD(5) * t374;
	t387 = t368 * t397;
	t386 = t369 * t395;
	t385 = t369 * t394;
	t378 = -t368 * t373 + t370 * t396;
	t355 = t378 * qJD(2);
	t382 = -t378 * t390 + t355;
	t361 = t368 * t396 + t370 * t373;
	t357 = t361 * qJD(2);
	t381 = t361 * t390 - t357;
	t380 = (qJD(2) - t390) * t375;
	t360 = t368 * t375 + t370 * t397;
	t349 = -t360 * t365 - t366 * t400;
	t379 = -t360 * t366 + t365 * t400;
	t362 = t370 * t375 - t387;
	t351 = -t362 * t365 + t366 * t401;
	t352 = t362 * t366 + t365 * t401;
	t354 = t371 * t365 + t366 * t399;
	t353 = -t365 * t399 + t371 * t366;
	t356 = t360 * qJD(2);
	t377 = qJD(5) * t360 - t356 * t366 - t378 * t393;
	t358 = -qJD(2) * t387 + t370 * t394;
	t376 = qJD(5) * t362 - t358 * t366 + t361 * t393;
	t348 = qJD(4) * t353 + t366 * t385;
	t347 = -qJD(4) * t354 - t365 * t385;
	t346 = qJD(4) * t351 - t357 * t366;
	t345 = -qJD(4) * t352 + t357 * t365;
	t344 = qJD(4) * t349 + t355 * t366;
	t343 = qJD(4) * t379 - t355 * t365;
	t1 = [0, t372 * t381 + t376 * t374, 0, t345 * t374 - t351 * t389, -t346 * t372 + t358 * t374 + (-t352 * t374 - t361 * t372) * qJD(5), 0; 0, t382 * t372 + t377 * t374, 0, t343 * t374 - t349 * t389, -t344 * t372 + t356 * t374 + (t372 * t378 + t374 * t379) * qJD(5), 0; 0, (t372 * t380 - t402 * t374) * t369, 0, t347 * t374 - t353 * t389, t374 * t386 - t348 * t372 + (-t354 * t374 + t372 * t398) * qJD(5), 0; 0, -t376 * t372 + t374 * t381, 0, -t345 * t372 - t351 * t388, -t346 * t374 - t358 * t372 + (t352 * t372 - t361 * t374) * qJD(5), 0; 0, -t377 * t372 + t382 * t374, 0, -t343 * t372 - t349 * t388, -t344 * t374 - t356 * t372 + (-t372 * t379 + t374 * t378) * qJD(5), 0; 0, (t402 * t372 + t374 * t380) * t369, 0, -t347 * t372 - t353 * t388, -t372 * t386 - t348 * t374 + (t354 * t372 + t374 * t398) * qJD(5), 0; 0, -t358 * t365 - t361 * t392, 0, t346, 0, 0; 0, -t356 * t365 + t378 * t392, 0, t344, 0, 0; 0, (-t365 * t395 + t366 * t391) * t369, 0, t348, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:21
	% EndTime: 2019-10-09 21:48:22
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (246->60), mult. (516->128), div. (0->0), fcn. (564->10), ass. (0->56)
	t447 = sin(pkin(10));
	t448 = sin(pkin(6));
	t480 = t447 * t448;
	t449 = cos(pkin(10));
	t479 = t448 * t449;
	t452 = sin(qJ(2));
	t478 = t448 * t452;
	t450 = cos(pkin(6));
	t477 = t450 * t452;
	t454 = cos(qJ(2));
	t476 = t450 * t454;
	t451 = sin(qJ(5));
	t475 = t451 * t454;
	t453 = cos(qJ(5));
	t474 = t453 * t454;
	t473 = qJD(2) * t452;
	t472 = qJD(2) * t454;
	t446 = pkin(11) + qJ(4);
	t444 = sin(t446);
	t471 = qJD(4) * t444;
	t445 = cos(t446);
	t470 = qJD(4) * t445;
	t469 = qJD(4) * t454;
	t468 = qJD(5) * t445;
	t467 = qJD(5) * t451;
	t466 = qJD(5) * t453;
	t465 = t447 * t477;
	t464 = t448 * t473;
	t463 = t448 * t472;
	t462 = -qJD(2) + t468;
	t458 = -t447 * t452 + t449 * t476;
	t434 = t458 * qJD(2);
	t461 = -t458 * t468 + t434;
	t440 = t447 * t476 + t449 * t452;
	t436 = t440 * qJD(2);
	t460 = t440 * t468 - t436;
	t439 = t447 * t454 + t449 * t477;
	t428 = -t439 * t444 - t445 * t479;
	t459 = -t439 * t445 + t444 * t479;
	t441 = t449 * t454 - t465;
	t430 = -t441 * t444 + t445 * t480;
	t431 = t441 * t445 + t444 * t480;
	t433 = t444 * t450 + t445 * t478;
	t432 = -t444 * t478 + t445 * t450;
	t435 = t439 * qJD(2);
	t457 = qJD(5) * t439 - t435 * t445 - t458 * t471;
	t437 = -qJD(2) * t465 + t449 * t472;
	t456 = qJD(5) * t441 - t437 * t445 + t440 * t471;
	t455 = -t444 * t469 + (-qJD(2) * t445 + qJD(5)) * t452;
	t427 = qJD(4) * t432 + t445 * t463;
	t426 = -qJD(4) * t433 - t444 * t463;
	t425 = qJD(4) * t430 - t436 * t445;
	t424 = -qJD(4) * t431 + t436 * t444;
	t423 = qJD(4) * t428 + t434 * t445;
	t422 = qJD(4) * t459 - t434 * t444;
	t1 = [0, t451 * t460 + t453 * t456, 0, t424 * t453 - t430 * t467, -t425 * t451 + t437 * t453 + (-t431 * t453 - t440 * t451) * qJD(5), 0; 0, t451 * t461 + t453 * t457, 0, t422 * t453 - t428 * t467, -t423 * t451 + t435 * t453 + (t451 * t458 + t453 * t459) * qJD(5), 0; 0, (t453 * t455 - t462 * t475) * t448, 0, t426 * t453 - t432 * t467, t453 * t464 - t427 * t451 + (-t433 * t453 + t448 * t475) * qJD(5), 0; 0, -t437 * t444 - t440 * t470, 0, t425, 0, 0; 0, -t435 * t444 + t458 * t470, 0, t423, 0, 0; 0, (-t444 * t473 + t445 * t469) * t448, 0, t427, 0, 0; 0, t451 * t456 - t453 * t460, 0, t424 * t451 + t430 * t466, t425 * t453 + t437 * t451 + (-t431 * t451 + t440 * t453) * qJD(5), 0; 0, t451 * t457 - t453 * t461, 0, t422 * t451 + t428 * t466, t423 * t453 + t435 * t451 + (t451 * t459 - t453 * t458) * qJD(5), 0; 0, (t451 * t455 + t462 * t474) * t448, 0, t426 * t451 + t432 * t466, t451 * t464 + t427 * t453 + (-t433 * t451 - t448 * t474) * qJD(5), 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end