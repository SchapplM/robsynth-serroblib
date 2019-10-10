% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
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
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(11));
	t17 = sin(pkin(11));
	t1 = [-t18 * t21, 0, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0, 0; t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:22
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t43 = sin(qJ(1));
	t48 = qJD(1) * t43;
	t44 = cos(qJ(1));
	t47 = qJD(1) * t44;
	t46 = qJD(3) * t43;
	t45 = qJD(3) * t44;
	t42 = pkin(11) + qJ(3);
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
	% StartTime: 2019-10-10 09:04:22
	% EndTime: 2019-10-10 09:04:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (89->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t72 = pkin(11) + qJ(3) + qJ(4);
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
	% StartTime: 2019-10-10 09:04:22
	% EndTime: 2019-10-10 09:04:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (181->17), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t91 = pkin(11) + qJ(3) + qJ(4) + qJ(5);
	t89 = sin(t91);
	t92 = qJD(3) + qJD(4) + qJD(5);
	t100 = t92 * t89;
	t90 = cos(t91);
	t99 = t92 * t90;
	t93 = sin(qJ(1));
	t98 = t92 * t93;
	t94 = cos(qJ(1));
	t97 = t92 * t94;
	t96 = qJD(1) * t93;
	t95 = qJD(1) * t94;
	t88 = t89 * t98 - t90 * t95;
	t87 = t89 * t95 + t90 * t98;
	t86 = t89 * t97 + t90 * t96;
	t85 = t89 * t96 - t90 * t97;
	t1 = [t88, 0, t85, t85, t85, 0; -t86, 0, -t87, -t87, -t87, 0; 0, 0, -t100, -t100, -t100, 0; t87, 0, t86, t86, t86, 0; t85, 0, t88, t88, t88, 0; 0, 0, -t99, -t99, -t99, 0; -t96, 0, 0, 0, 0, 0; t95, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:23
	% EndTime: 2019-10-10 09:04:24
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (435->29), mult. (279->57), div. (0->0), fcn. (279->6), ass. (0->41)
	t352 = qJD(3) + qJD(4) + qJD(5);
	t353 = sin(qJ(6));
	t375 = t352 * t353;
	t354 = sin(qJ(1));
	t374 = t352 * t354;
	t355 = cos(qJ(6));
	t373 = t352 * t355;
	t356 = cos(qJ(1));
	t372 = t352 * t356;
	t371 = t355 * t356;
	t370 = qJD(1) * t354;
	t369 = qJD(1) * t356;
	t368 = qJD(6) * t353;
	t367 = qJD(6) * t355;
	t366 = qJD(6) * t356;
	t351 = pkin(11) + qJ(3) + qJ(4) + qJ(5);
	t349 = sin(t351);
	t365 = t349 * t373;
	t364 = t349 * t374;
	t350 = cos(t351);
	t363 = t350 * t374;
	t362 = t349 * t372;
	t361 = t350 * t372;
	t360 = qJD(6) * t350 - qJD(1);
	t359 = qJD(1) * t350 - qJD(6);
	t358 = t360 * t353;
	t357 = t359 * t354 + t362;
	t348 = t352 * t350;
	t347 = t350 * t369 - t364;
	t346 = -t350 * t370 - t362;
	t345 = -t350 * t368 - t365;
	t344 = t349 * t375 - t350 * t367;
	t343 = -t355 * t363 + (t354 * t368 - t355 * t369) * t349;
	t342 = t353 * t363 + (t353 * t369 + t354 * t367) * t349;
	t341 = -t355 * t361 + (t353 * t366 + t355 * t370) * t349;
	t340 = t353 * t361 + (-t353 * t370 + t355 * t366) * t349;
	t339 = -t359 * t371 + (t358 + t365) * t354;
	t338 = t360 * t355 * t354 + (t359 * t356 - t364) * t353;
	t337 = t357 * t355 + t356 * t358;
	t336 = t357 * t353 - t360 * t371;
	t1 = [t339, 0, t341, t341, t341, t336; -t337, 0, t343, t343, t343, -t338; 0, 0, t345, t345, t345, -t349 * t367 - t350 * t375; t338, 0, t340, t340, t340, t337; t336, 0, t342, t342, t342, t339; 0, 0, t344, t344, t344, t349 * t368 - t350 * t373; -t349 * t369 - t363, 0, t346, t346, t346, 0; -t349 * t370 + t361, 0, t347, t347, t347, 0; 0, 0, t348, t348, t348, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end