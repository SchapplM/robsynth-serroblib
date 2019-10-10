% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
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
	% StartTime: 2019-10-09 22:14:36
	% EndTime: 2019-10-09 22:14:36
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
	t210 = sin(pkin(10));
	t212 = cos(pkin(10));
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
	% StartTime: 2019-10-09 22:14:36
	% EndTime: 2019-10-09 22:14:36
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t257 = sin(pkin(6));
	t260 = sin(qJ(3));
	t273 = t257 * t260;
	t262 = cos(qJ(3));
	t272 = t257 * t262;
	t259 = cos(pkin(6));
	t261 = sin(qJ(2));
	t271 = t259 * t261;
	t263 = cos(qJ(2));
	t270 = t259 * t263;
	t269 = qJD(2) * t261;
	t268 = qJD(3) * t260;
	t267 = qJD(3) * t262;
	t266 = qJD(3) * t263;
	t265 = t257 * qJD(2) * t263;
	t256 = sin(pkin(10));
	t258 = cos(pkin(10));
	t252 = -t256 * t261 + t258 * t270;
	t253 = t256 * t263 + t258 * t271;
	t254 = -t256 * t270 - t258 * t261;
	t264 = t256 * t271 - t258 * t263;
	t251 = t264 * qJD(2);
	t250 = t254 * qJD(2);
	t249 = t253 * qJD(2);
	t248 = t252 * qJD(2);
	t1 = [0, t250, 0, 0, 0, 0; 0, t248, 0, 0, 0, 0; 0, t265, 0, 0, 0, 0; 0, -t251 * t262 + t254 * t268, t250 * t260 + (t256 * t273 - t262 * t264) * qJD(3), 0, 0, 0; 0, t249 * t262 + t252 * t268, t248 * t260 + (t253 * t262 - t258 * t273) * qJD(3), 0, 0, 0; 0, (t260 * t266 + t262 * t269) * t257, t260 * t265 + (t259 * t260 + t261 * t272) * qJD(3), 0, 0, 0; 0, t251 * t260 + t254 * t267, t250 * t262 + (t256 * t272 + t260 * t264) * qJD(3), 0, 0, 0; 0, -t249 * t260 + t252 * t267, t248 * t262 + (-t253 * t260 - t258 * t272) * qJD(3), 0, 0, 0; 0, (-t260 * t269 + t262 * t266) * t257, t262 * t265 + (t259 * t262 - t261 * t273) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:37
	% EndTime: 2019-10-09 22:14:37
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (61->29), mult. (240->82), div. (0->0), fcn. (252->10), ass. (0->35)
	t306 = sin(pkin(6));
	t310 = sin(qJ(3));
	t327 = t306 * t310;
	t312 = cos(qJ(3));
	t326 = t306 * t312;
	t309 = cos(pkin(6));
	t311 = sin(qJ(2));
	t325 = t309 * t311;
	t313 = cos(qJ(2));
	t324 = t309 * t313;
	t323 = t310 * t311;
	t322 = t311 * t312;
	t321 = qJD(3) * t310;
	t320 = qJD(3) * t312;
	t319 = qJD(3) * t313;
	t318 = qJD(2) * t306 * t313;
	t317 = t312 * t319;
	t305 = sin(pkin(10));
	t308 = cos(pkin(10));
	t300 = -t305 * t311 + t308 * t324;
	t301 = t305 * t313 + t308 * t325;
	t302 = -t305 * t324 - t308 * t311;
	t316 = t305 * t325 - t308 * t313;
	t297 = t301 * qJD(2);
	t315 = -t297 * t310 + t300 * t320;
	t299 = t316 * qJD(2);
	t314 = t299 * t310 + t302 * t320;
	t307 = cos(pkin(11));
	t304 = sin(pkin(11));
	t298 = t302 * qJD(2);
	t296 = t300 * qJD(2);
	t295 = t312 * t318 + (-t306 * t323 + t309 * t312) * qJD(3);
	t294 = t298 * t312 + (t305 * t326 + t310 * t316) * qJD(3);
	t293 = t296 * t312 + (-t301 * t310 - t308 * t326) * qJD(3);
	t1 = [0, t298 * t307 + t314 * t304, t294 * t304, 0, 0, 0; 0, t296 * t307 + t315 * t304, t293 * t304, 0, 0, 0; 0, (t304 * t317 + (-t304 * t323 + t307 * t313) * qJD(2)) * t306, t295 * t304, 0, 0, 0; 0, -t298 * t304 + t314 * t307, t294 * t307, 0, 0, 0; 0, -t296 * t304 + t315 * t307, t293 * t307, 0, 0, 0; 0, (t307 * t317 + (-t304 * t313 - t307 * t323) * qJD(2)) * t306, t295 * t307, 0, 0, 0; 0, t299 * t312 - t302 * t321, -t298 * t310 + (-t305 * t327 + t312 * t316) * qJD(3), 0, 0, 0; 0, -t297 * t312 - t300 * t321, -t296 * t310 + (-t301 * t312 + t308 * t327) * qJD(3), 0, 0, 0; 0, (-qJD(2) * t322 - t310 * t319) * t306, -t310 * t318 + (-t306 * t322 - t309 * t310) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:37
	% EndTime: 2019-10-09 22:14:37
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (219->63), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->54)
	t381 = cos(qJ(2));
	t378 = sin(qJ(3));
	t393 = qJD(6) * t378;
	t406 = t381 * (qJD(2) + t393);
	t375 = sin(pkin(6));
	t405 = t375 * t378;
	t380 = cos(qJ(3));
	t404 = t375 * t380;
	t403 = t375 * t381;
	t377 = cos(pkin(6));
	t379 = sin(qJ(2));
	t402 = t377 * t379;
	t401 = t377 * t381;
	t400 = qJD(2) * t379;
	t399 = qJD(2) * t381;
	t398 = qJD(3) * t378;
	t397 = qJD(3) * t380;
	t396 = qJD(3) * t381;
	t373 = pkin(11) + qJ(6);
	t371 = sin(t373);
	t395 = qJD(6) * t371;
	t372 = cos(t373);
	t394 = qJD(6) * t372;
	t374 = sin(pkin(10));
	t392 = t374 * t402;
	t391 = t375 * t400;
	t390 = t375 * t399;
	t376 = cos(pkin(10));
	t385 = -t374 * t379 + t376 * t401;
	t359 = t385 * qJD(2);
	t388 = -t385 * t393 - t359;
	t365 = t374 * t401 + t376 * t379;
	t361 = t365 * qJD(2);
	t387 = t365 * t393 + t361;
	t364 = t374 * t381 + t376 * t402;
	t353 = t364 * t378 + t376 * t404;
	t354 = t364 * t380 - t376 * t405;
	t366 = t376 * t381 - t392;
	t386 = -t366 * t378 + t374 * t404;
	t356 = t366 * t380 + t374 * t405;
	t368 = t377 * t378 + t379 * t404;
	t367 = -t377 * t380 + t379 * t405;
	t360 = t364 * qJD(2);
	t384 = -qJD(6) * t364 - t360 * t378 + t385 * t397;
	t362 = -qJD(2) * t392 + t376 * t399;
	t383 = -qJD(6) * t366 - t362 * t378 - t365 * t397;
	t382 = t380 * t396 + (-qJD(2) * t378 - qJD(6)) * t379;
	t358 = -t367 * qJD(3) + t380 * t390;
	t357 = t368 * qJD(3) + t378 * t390;
	t352 = t386 * qJD(3) - t361 * t380;
	t351 = t356 * qJD(3) - t361 * t378;
	t350 = -t353 * qJD(3) + t359 * t380;
	t349 = t354 * qJD(3) + t359 * t378;
	t1 = [0, t383 * t371 - t387 * t372, t352 * t371 + t356 * t394, 0, 0, t351 * t372 - t362 * t371 + (-t365 * t372 + t371 * t386) * qJD(6); 0, t384 * t371 - t388 * t372, t350 * t371 + t354 * t394, 0, 0, t349 * t372 - t360 * t371 + (-t353 * t371 + t372 * t385) * qJD(6); 0, (t382 * t371 + t372 * t406) * t375, t358 * t371 + t368 * t394, 0, 0, -t371 * t391 + t357 * t372 + (-t367 * t371 + t372 * t403) * qJD(6); 0, t387 * t371 + t383 * t372, t352 * t372 - t356 * t395, 0, 0, -t351 * t371 - t362 * t372 + (t365 * t371 + t372 * t386) * qJD(6); 0, t388 * t371 + t384 * t372, t350 * t372 - t354 * t395, 0, 0, -t349 * t371 - t360 * t372 + (-t353 * t372 - t371 * t385) * qJD(6); 0, (-t371 * t406 + t382 * t372) * t375, t358 * t372 - t368 * t395, 0, 0, -t372 * t391 - t357 * t371 + (-t367 * t372 - t371 * t403) * qJD(6); 0, -t362 * t380 + t365 * t398, -t351, 0, 0, 0; 0, -t360 * t380 - t385 * t398, -t349, 0, 0, 0; 0, (-t378 * t396 - t380 * t400) * t375, -t357, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end