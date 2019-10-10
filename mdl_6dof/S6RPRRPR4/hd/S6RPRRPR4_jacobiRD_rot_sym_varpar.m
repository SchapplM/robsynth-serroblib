% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
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
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(10));
	t17 = sin(pkin(10));
	t1 = [-t18 * t21, 0, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0, 0; t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t43 = sin(qJ(1));
	t48 = qJD(1) * t43;
	t44 = cos(qJ(1));
	t47 = qJD(1) * t44;
	t46 = qJD(3) * t43;
	t45 = qJD(3) * t44;
	t42 = pkin(10) + qJ(3);
	t41 = cos(t42);
	t40 = sin(t42);
	t39 = t40 * t46 - t41 * t47;
	t38 = t40 * t47 + t41 * t46;
	t37 = t40 * t45 + t41 * t48;
	t36 = t40 * t48 - t41 * t45;
	t1 = [t39, 0, t36, 0, 0, 0; -t37, 0, -t38, 0, 0, 0; 0, 0, -qJD(3) * t40, 0, 0, 0; t38, 0, t37, 0, 0, 0; t36, 0, t39, 0, 0, 0; 0, 0, -qJD(3) * t41, 0, 0, 0; -t48, 0, 0, 0, 0, 0; t47, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (89->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t72 = pkin(10) + qJ(3) + qJ(4);
	t70 = sin(t72);
	t73 = qJD(3) + qJD(4);
	t81 = t73 * t70;
	t71 = cos(t72);
	t80 = t73 * t71;
	t74 = sin(qJ(1));
	t79 = t73 * t74;
	t75 = cos(qJ(1));
	t78 = t73 * t75;
	t77 = qJD(1) * t74;
	t76 = qJD(1) * t75;
	t69 = t70 * t79 - t71 * t76;
	t68 = t70 * t76 + t71 * t79;
	t67 = t70 * t78 + t71 * t77;
	t66 = t70 * t77 - t71 * t78;
	t1 = [t69, 0, t66, t66, 0, 0; -t67, 0, -t68, -t68, 0, 0; 0, 0, -t81, -t81, 0, 0; t68, 0, t67, t67, 0, 0; t66, 0, t69, t69, 0, 0; 0, 0, -t80, -t80, 0, 0; -t77, 0, 0, 0, 0, 0; t76, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:47
	% EndTime: 2019-10-10 01:28:48
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (132->22), mult. (114->36), div. (0->0), fcn. (114->6), ass. (0->31)
	t283 = pkin(10) + qJ(3) + qJ(4);
	t281 = sin(t283);
	t284 = qJD(3) + qJD(4);
	t302 = t281 * t284;
	t287 = sin(qJ(1));
	t301 = t284 * t287;
	t288 = cos(qJ(1));
	t300 = t284 * t288;
	t285 = sin(pkin(11));
	t299 = t285 * t287;
	t298 = t285 * t288;
	t286 = cos(pkin(11));
	t297 = t286 * t287;
	t296 = t286 * t288;
	t295 = qJD(1) * t287;
	t294 = qJD(1) * t288;
	t293 = t286 * t302;
	t292 = t281 * t301;
	t291 = t281 * t300;
	t282 = cos(t283);
	t290 = t281 * t294 + t282 * t301;
	t289 = t281 * t295 - t282 * t300;
	t280 = t284 * t282;
	t279 = t285 * t302;
	t278 = t282 * t294 - t292;
	t277 = -t282 * t295 - t291;
	t276 = t290 * t286;
	t275 = t290 * t285;
	t274 = t289 * t286;
	t273 = t289 * t285;
	t1 = [t286 * t292 + (-t282 * t296 - t299) * qJD(1), 0, t274, t274, 0, 0; -t286 * t291 + (-t282 * t297 + t298) * qJD(1), 0, -t276, -t276, 0, 0; 0, 0, -t293, -t293, 0, 0; -t285 * t292 + (t282 * t298 - t297) * qJD(1), 0, -t273, -t273, 0, 0; t285 * t291 + (t282 * t299 + t296) * qJD(1), 0, t275, t275, 0, 0; 0, 0, t279, t279, 0, 0; -t290, 0, t277, t277, 0, 0; -t289, 0, t278, t278, 0, 0; 0, 0, t280, t280, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:48
	% EndTime: 2019-10-10 01:28:48
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (314->29), mult. (226->56), div. (0->0), fcn. (226->6), ass. (0->42)
	t349 = cos(qJ(1));
	t345 = pkin(10) + qJ(3) + qJ(4);
	t342 = cos(t345);
	t353 = qJD(6) * t342 - qJD(1);
	t368 = t349 * t353;
	t352 = qJD(1) * t342 - qJD(6);
	t341 = sin(t345);
	t347 = qJD(3) + qJD(4);
	t348 = sin(qJ(1));
	t365 = t347 * t348;
	t357 = t341 * t365;
	t367 = t352 * t349 - t357;
	t366 = t341 * t347;
	t340 = t347 * t342;
	t364 = t347 * t349;
	t363 = qJD(1) * t348;
	t362 = qJD(1) * t349;
	t346 = pkin(11) + qJ(6);
	t343 = sin(t346);
	t361 = qJD(6) * t343;
	t344 = cos(t346);
	t360 = qJD(6) * t344;
	t359 = qJD(6) * t348;
	t358 = qJD(6) * t349;
	t356 = t342 * t365;
	t355 = t341 * t364;
	t354 = t342 * t364;
	t351 = t353 * t348;
	t350 = t352 * t348 + t355;
	t339 = t342 * t362 - t357;
	t338 = -t342 * t363 - t355;
	t337 = -t342 * t361 - t344 * t366;
	t336 = -t342 * t360 + t343 * t366;
	t335 = -t344 * t356 + (t343 * t359 - t344 * t362) * t341;
	t334 = t343 * t356 + (t343 * t362 + t344 * t359) * t341;
	t333 = -t344 * t354 + (t343 * t358 + t344 * t363) * t341;
	t332 = t343 * t354 + (-t343 * t363 + t344 * t358) * t341;
	t331 = t343 * t351 - t367 * t344;
	t330 = t367 * t343 + t344 * t351;
	t329 = t343 * t368 + t350 * t344;
	t328 = t350 * t343 - t344 * t368;
	t1 = [t331, 0, t333, t333, 0, t328; -t329, 0, t335, t335, 0, -t330; 0, 0, t337, t337, 0, -t343 * t340 - t341 * t360; t330, 0, t332, t332, 0, t329; t328, 0, t334, t334, 0, t331; 0, 0, t336, t336, 0, -t344 * t340 + t341 * t361; -t341 * t362 - t356, 0, t338, t338, 0, 0; -t341 * t363 + t354, 0, t339, t339, 0, 0; 0, 0, t340, t340, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end