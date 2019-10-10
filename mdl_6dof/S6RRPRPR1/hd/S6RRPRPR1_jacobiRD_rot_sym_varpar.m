% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
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
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t50 = sin(qJ(1));
	t55 = qJD(1) * t50;
	t51 = cos(qJ(1));
	t54 = qJD(1) * t51;
	t53 = qJD(2) * t50;
	t52 = qJD(2) * t51;
	t49 = qJ(2) + pkin(10);
	t48 = cos(t49);
	t47 = sin(t49);
	t46 = t47 * t53 - t48 * t54;
	t45 = t47 * t54 + t48 * t53;
	t44 = t47 * t52 + t48 * t55;
	t43 = t47 * t55 - t48 * t52;
	t1 = [t46, t43, 0, 0, 0, 0; -t44, -t45, 0, 0, 0, 0; 0, -qJD(2) * t47, 0, 0, 0, 0; t45, t44, 0, 0, 0, 0; t43, t46, 0, 0, 0, 0; 0, -qJD(2) * t48, 0, 0, 0, 0; -t55, 0, 0, 0, 0, 0; t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (89->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t77 = qJ(2) + pkin(10) + qJ(4);
	t75 = sin(t77);
	t78 = qJD(2) + qJD(4);
	t86 = t78 * t75;
	t76 = cos(t77);
	t85 = t78 * t76;
	t79 = sin(qJ(1));
	t84 = t78 * t79;
	t80 = cos(qJ(1));
	t83 = t78 * t80;
	t82 = qJD(1) * t79;
	t81 = qJD(1) * t80;
	t74 = t75 * t84 - t76 * t81;
	t73 = t75 * t81 + t76 * t84;
	t72 = t75 * t83 + t76 * t82;
	t71 = t75 * t82 - t76 * t83;
	t1 = [t74, t71, 0, t71, 0, 0; -t72, -t73, 0, -t73, 0, 0; 0, -t86, 0, -t86, 0, 0; t73, t72, 0, t72, 0, 0; t71, t74, 0, t74, 0, 0; 0, -t85, 0, -t85, 0, 0; -t82, 0, 0, 0, 0, 0; t81, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:27
	% EndTime: 2019-10-10 10:04:27
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (132->22), mult. (114->36), div. (0->0), fcn. (114->6), ass. (0->31)
	t288 = qJ(2) + pkin(10) + qJ(4);
	t286 = sin(t288);
	t289 = qJD(2) + qJD(4);
	t307 = t286 * t289;
	t292 = sin(qJ(1));
	t306 = t289 * t292;
	t293 = cos(qJ(1));
	t305 = t289 * t293;
	t290 = sin(pkin(11));
	t304 = t290 * t292;
	t303 = t290 * t293;
	t291 = cos(pkin(11));
	t302 = t291 * t292;
	t301 = t291 * t293;
	t300 = qJD(1) * t292;
	t299 = qJD(1) * t293;
	t298 = t291 * t307;
	t297 = t286 * t306;
	t296 = t286 * t305;
	t287 = cos(t288);
	t295 = t286 * t299 + t287 * t306;
	t294 = t286 * t300 - t287 * t305;
	t285 = t289 * t287;
	t284 = t290 * t307;
	t283 = t287 * t299 - t297;
	t282 = -t287 * t300 - t296;
	t281 = t295 * t291;
	t280 = t295 * t290;
	t279 = t294 * t291;
	t278 = t294 * t290;
	t1 = [t291 * t297 + (-t287 * t301 - t304) * qJD(1), t279, 0, t279, 0, 0; -t291 * t296 + (-t287 * t302 + t303) * qJD(1), -t281, 0, -t281, 0, 0; 0, -t298, 0, -t298, 0, 0; -t290 * t297 + (t287 * t303 - t302) * qJD(1), -t278, 0, -t278, 0, 0; t290 * t296 + (t287 * t304 + t301) * qJD(1), t280, 0, t280, 0, 0; 0, t284, 0, t284, 0, 0; -t295, t282, 0, t282, 0, 0; -t294, t283, 0, t283, 0, 0; 0, t285, 0, t285, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:27
	% EndTime: 2019-10-10 10:04:28
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (314->29), mult. (226->56), div. (0->0), fcn. (226->6), ass. (0->42)
	t354 = cos(qJ(1));
	t350 = qJ(2) + pkin(10) + qJ(4);
	t347 = cos(t350);
	t358 = qJD(6) * t347 - qJD(1);
	t373 = t354 * t358;
	t357 = qJD(1) * t347 - qJD(6);
	t346 = sin(t350);
	t352 = qJD(2) + qJD(4);
	t353 = sin(qJ(1));
	t370 = t352 * t353;
	t362 = t346 * t370;
	t372 = t357 * t354 - t362;
	t371 = t346 * t352;
	t345 = t352 * t347;
	t369 = t352 * t354;
	t368 = qJD(1) * t353;
	t367 = qJD(1) * t354;
	t351 = pkin(11) + qJ(6);
	t348 = sin(t351);
	t366 = qJD(6) * t348;
	t349 = cos(t351);
	t365 = qJD(6) * t349;
	t364 = qJD(6) * t353;
	t363 = qJD(6) * t354;
	t361 = t347 * t370;
	t360 = t346 * t369;
	t359 = t347 * t369;
	t356 = t358 * t353;
	t355 = t357 * t353 + t360;
	t344 = t347 * t367 - t362;
	t343 = -t347 * t368 - t360;
	t342 = -t347 * t366 - t349 * t371;
	t341 = -t347 * t365 + t348 * t371;
	t340 = -t349 * t361 + (t348 * t364 - t349 * t367) * t346;
	t339 = t348 * t361 + (t348 * t367 + t349 * t364) * t346;
	t338 = -t349 * t359 + (t348 * t363 + t349 * t368) * t346;
	t337 = t348 * t359 + (-t348 * t368 + t349 * t363) * t346;
	t336 = t348 * t356 - t372 * t349;
	t335 = t372 * t348 + t349 * t356;
	t334 = t348 * t373 + t355 * t349;
	t333 = t355 * t348 - t349 * t373;
	t1 = [t336, t338, 0, t338, 0, t333; -t334, t340, 0, t340, 0, -t335; 0, t342, 0, t342, 0, -t348 * t345 - t346 * t365; t335, t337, 0, t337, 0, t334; t333, t339, 0, t339, 0, t336; 0, t341, 0, t341, 0, -t349 * t345 + t346 * t366; -t346 * t367 - t361, t343, 0, t343, 0, 0; -t346 * t368 + t359, t344, 0, t344, 0, 0; 0, t345, 0, t345, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end