% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRP3
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
% Datum: 2019-10-09 21:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRRP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:46:28
	% EndTime: 2019-10-09 21:46:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:46:28
	% EndTime: 2019-10-09 21:46:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:46:28
	% EndTime: 2019-10-09 21:46:28
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
	% StartTime: 2019-10-09 21:46:29
	% EndTime: 2019-10-09 21:46:29
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
	% StartTime: 2019-10-09 21:46:29
	% EndTime: 2019-10-09 21:46:29
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-09 21:46:30
	% EndTime: 2019-10-09 21:46:31
	% DurationCPUTime: 0.38s
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
	t1 = [0, t381 * t372 + t376 * t374, 0, t345 * t374 - t351 * t389, -t346 * t372 + t358 * t374 + (-t352 * t374 - t361 * t372) * qJD(5), 0; 0, t382 * t372 + t377 * t374, 0, t343 * t374 - t349 * t389, -t344 * t372 + t356 * t374 + (t372 * t378 + t374 * t379) * qJD(5), 0; 0, (t372 * t380 - t402 * t374) * t369, 0, t347 * t374 - t353 * t389, t374 * t386 - t348 * t372 + (-t354 * t374 + t372 * t398) * qJD(5), 0; 0, -t376 * t372 + t381 * t374, 0, -t345 * t372 - t351 * t388, -t346 * t374 - t358 * t372 + (t352 * t372 - t361 * t374) * qJD(5), 0; 0, -t377 * t372 + t382 * t374, 0, -t343 * t372 - t349 * t388, -t344 * t374 - t356 * t372 + (-t372 * t379 + t374 * t378) * qJD(5), 0; 0, (t402 * t372 + t374 * t380) * t369, 0, -t347 * t372 - t353 * t388, -t372 * t386 - t348 * t374 + (t354 * t372 + t374 * t398) * qJD(5), 0; 0, -t358 * t365 - t361 * t392, 0, t346, 0, 0; 0, -t356 * t365 + t378 * t392, 0, t344, 0, 0; 0, (-t365 * t395 + t366 * t391) * t369, 0, t348, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:46:31
	% EndTime: 2019-10-09 21:46:31
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (246->60), mult. (516->128), div. (0->0), fcn. (564->10), ass. (0->55)
	t390 = pkin(11) + qJ(4);
	t388 = sin(t390);
	t389 = cos(t390);
	t396 = sin(qJ(2));
	t398 = cos(qJ(2));
	t414 = qJD(4) * t398;
	t425 = (qJD(2) * t389 - qJD(5)) * t396 + t388 * t414;
	t391 = sin(pkin(10));
	t392 = sin(pkin(6));
	t424 = t391 * t392;
	t393 = cos(pkin(10));
	t423 = t392 * t393;
	t422 = t392 * t396;
	t421 = t392 * t398;
	t394 = cos(pkin(6));
	t420 = t394 * t396;
	t419 = t394 * t398;
	t418 = qJD(2) * t396;
	t417 = qJD(2) * t398;
	t416 = qJD(4) * t388;
	t415 = qJD(4) * t389;
	t413 = qJD(5) * t389;
	t395 = sin(qJ(5));
	t412 = qJD(5) * t395;
	t397 = cos(qJ(5));
	t411 = qJD(5) * t397;
	t410 = t391 * t420;
	t409 = t392 * t418;
	t408 = t392 * t417;
	t401 = -t391 * t396 + t393 * t419;
	t378 = t401 * qJD(2);
	t405 = -t401 * t413 + t378;
	t384 = t391 * t419 + t393 * t396;
	t380 = t384 * qJD(2);
	t404 = t384 * t413 - t380;
	t403 = (qJD(2) - t413) * t398;
	t383 = t391 * t398 + t393 * t420;
	t372 = -t383 * t388 - t389 * t423;
	t402 = -t383 * t389 + t388 * t423;
	t385 = t393 * t398 - t410;
	t374 = -t385 * t388 + t389 * t424;
	t375 = t385 * t389 + t388 * t424;
	t377 = t394 * t388 + t389 * t422;
	t376 = -t388 * t422 + t394 * t389;
	t379 = t383 * qJD(2);
	t400 = qJD(5) * t383 - t379 * t389 - t401 * t416;
	t381 = -qJD(2) * t410 + t393 * t417;
	t399 = qJD(5) * t385 - t381 * t389 + t384 * t416;
	t371 = t376 * qJD(4) + t389 * t408;
	t370 = -t377 * qJD(4) - t388 * t408;
	t369 = t374 * qJD(4) - t380 * t389;
	t368 = -t375 * qJD(4) + t380 * t388;
	t367 = t372 * qJD(4) + t378 * t389;
	t366 = t402 * qJD(4) - t378 * t388;
	t1 = [0, t404 * t395 + t399 * t397, 0, t368 * t397 - t374 * t412, -t369 * t395 + t381 * t397 + (-t375 * t397 - t384 * t395) * qJD(5), 0; 0, t405 * t395 + t400 * t397, 0, t366 * t397 - t372 * t412, -t367 * t395 + t379 * t397 + (t395 * t401 + t397 * t402) * qJD(5), 0; 0, (t395 * t403 - t425 * t397) * t392, 0, t370 * t397 - t376 * t412, t397 * t409 - t371 * t395 + (-t377 * t397 + t395 * t421) * qJD(5), 0; 0, -t399 * t395 + t404 * t397, 0, -t368 * t395 - t374 * t411, -t369 * t397 - t381 * t395 + (t375 * t395 - t384 * t397) * qJD(5), 0; 0, -t400 * t395 + t405 * t397, 0, -t366 * t395 - t372 * t411, -t367 * t397 - t379 * t395 + (-t395 * t402 + t397 * t401) * qJD(5), 0; 0, (t425 * t395 + t397 * t403) * t392, 0, -t370 * t395 - t376 * t411, -t395 * t409 - t371 * t397 + (t377 * t395 + t397 * t421) * qJD(5), 0; 0, -t381 * t388 - t384 * t415, 0, t369, 0, 0; 0, -t379 * t388 + t401 * t415, 0, t367, 0, 0; 0, (-t388 * t418 + t389 * t414) * t392, 0, t371, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end