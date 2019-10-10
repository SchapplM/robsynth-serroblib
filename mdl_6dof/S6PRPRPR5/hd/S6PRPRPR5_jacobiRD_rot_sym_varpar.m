% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
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
	% StartTime: 2019-10-09 21:37:22
	% EndTime: 2019-10-09 21:37:22
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
	% StartTime: 2019-10-09 21:37:22
	% EndTime: 2019-10-09 21:37:22
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
	% StartTime: 2019-10-09 21:37:23
	% EndTime: 2019-10-09 21:37:23
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (66->23), mult. (140->60), div. (0->0), fcn. (148->8), ass. (0->28)
	t273 = sin(pkin(10));
	t274 = sin(pkin(6));
	t289 = t273 * t274;
	t275 = cos(pkin(10));
	t288 = t274 * t275;
	t277 = sin(qJ(2));
	t287 = t274 * t277;
	t276 = cos(pkin(6));
	t286 = t276 * t277;
	t278 = cos(qJ(2));
	t285 = t276 * t278;
	t284 = qJD(2) * t277;
	t272 = pkin(11) + qJ(4);
	t270 = sin(t272);
	t283 = qJD(4) * t270;
	t271 = cos(t272);
	t282 = qJD(4) * t271;
	t281 = qJD(4) * t278;
	t280 = t274 * qJD(2) * t278;
	t266 = -t273 * t277 + t275 * t285;
	t267 = t273 * t278 + t275 * t286;
	t268 = -t273 * t285 - t275 * t277;
	t279 = t273 * t286 - t275 * t278;
	t265 = t279 * qJD(2);
	t264 = t268 * qJD(2);
	t263 = t267 * qJD(2);
	t262 = t266 * qJD(2);
	t1 = [0, t264, 0, 0, 0, 0; 0, t262, 0, 0, 0, 0; 0, t280, 0, 0, 0, 0; 0, -t265 * t271 + t268 * t283, 0, t264 * t270 + (t270 * t289 - t271 * t279) * qJD(4), 0, 0; 0, t263 * t271 + t266 * t283, 0, t262 * t270 + (t267 * t271 - t270 * t288) * qJD(4), 0, 0; 0, (t270 * t281 + t271 * t284) * t274, 0, t270 * t280 + (t270 * t276 + t271 * t287) * qJD(4), 0, 0; 0, t265 * t270 + t268 * t282, 0, t264 * t271 + (t270 * t279 + t271 * t289) * qJD(4), 0, 0; 0, -t263 * t270 + t266 * t282, 0, t262 * t271 + (-t267 * t270 - t271 * t288) * qJD(4), 0, 0; 0, (-t270 * t284 + t271 * t281) * t274, 0, t271 * t280 + (-t270 * t287 + t271 * t276) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:24
	% EndTime: 2019-10-09 21:37:24
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (246->63), mult. (516->128), div. (0->0), fcn. (564->10), ass. (0->56)
	t376 = sin(pkin(10));
	t377 = sin(pkin(6));
	t409 = t376 * t377;
	t378 = cos(pkin(10));
	t408 = t377 * t378;
	t381 = sin(qJ(2));
	t407 = t377 * t381;
	t379 = cos(pkin(6));
	t406 = t379 * t381;
	t383 = cos(qJ(2));
	t405 = t379 * t383;
	t380 = sin(qJ(6));
	t404 = t380 * t383;
	t382 = cos(qJ(6));
	t403 = t382 * t383;
	t402 = qJD(2) * t381;
	t401 = qJD(2) * t383;
	t375 = pkin(11) + qJ(4);
	t373 = sin(t375);
	t400 = qJD(4) * t373;
	t374 = cos(t375);
	t399 = qJD(4) * t374;
	t398 = qJD(4) * t383;
	t397 = qJD(6) * t373;
	t396 = qJD(6) * t380;
	t395 = qJD(6) * t382;
	t394 = t376 * t406;
	t393 = t377 * t402;
	t392 = t377 * t401;
	t391 = qJD(2) + t397;
	t387 = -t376 * t381 + t378 * t405;
	t363 = t387 * qJD(2);
	t390 = -t387 * t397 - t363;
	t369 = t376 * t405 + t378 * t381;
	t365 = t369 * qJD(2);
	t389 = t369 * t397 + t365;
	t368 = t376 * t383 + t378 * t406;
	t357 = t368 * t373 + t374 * t408;
	t358 = t368 * t374 - t373 * t408;
	t370 = t378 * t383 - t394;
	t388 = -t370 * t373 + t374 * t409;
	t360 = t370 * t374 + t373 * t409;
	t362 = t379 * t373 + t374 * t407;
	t361 = t373 * t407 - t379 * t374;
	t364 = t368 * qJD(2);
	t386 = -qJD(6) * t368 - t364 * t373 + t387 * t399;
	t366 = -qJD(2) * t394 + t378 * t401;
	t385 = -qJD(6) * t370 - t366 * t373 - t369 * t399;
	t384 = t374 * t398 + (-qJD(2) * t373 - qJD(6)) * t381;
	t356 = -t361 * qJD(4) + t374 * t392;
	t355 = t362 * qJD(4) + t373 * t392;
	t354 = t388 * qJD(4) - t365 * t374;
	t353 = t360 * qJD(4) - t365 * t373;
	t352 = -t357 * qJD(4) + t363 * t374;
	t351 = t358 * qJD(4) + t363 * t373;
	t1 = [0, t385 * t380 - t389 * t382, 0, t354 * t380 + t360 * t395, 0, t353 * t382 - t366 * t380 + (-t369 * t382 + t380 * t388) * qJD(6); 0, t386 * t380 - t390 * t382, 0, t352 * t380 + t358 * t395, 0, t351 * t382 - t364 * t380 + (-t357 * t380 + t382 * t387) * qJD(6); 0, (t384 * t380 + t391 * t403) * t377, 0, t356 * t380 + t362 * t395, 0, -t380 * t393 + t355 * t382 + (-t361 * t380 + t377 * t403) * qJD(6); 0, t389 * t380 + t385 * t382, 0, t354 * t382 - t360 * t396, 0, -t353 * t380 - t366 * t382 + (t369 * t380 + t382 * t388) * qJD(6); 0, t390 * t380 + t386 * t382, 0, t352 * t382 - t358 * t396, 0, -t351 * t380 - t364 * t382 + (-t357 * t382 - t380 * t387) * qJD(6); 0, (t384 * t382 - t391 * t404) * t377, 0, t356 * t382 - t362 * t396, 0, -t382 * t393 - t355 * t380 + (-t361 * t382 - t377 * t404) * qJD(6); 0, -t366 * t374 + t369 * t400, 0, -t353, 0, 0; 0, -t364 * t374 - t387 * t400, 0, -t351, 0, 0; 0, (-t373 * t398 - t374 * t402) * t377, 0, -t355, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end