% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
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
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(10);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(10);
	t39 = cos(t40);
	t38 = sin(t40);
	t37 = t38 * t44 - t39 * t45;
	t36 = t38 * t43 + t39 * t46;
	t35 = t38 * t45 + t39 * t44;
	t34 = t38 * t46 - t39 * t43;
	t1 = [t37, 0, t34, 0, 0, 0; -t35, 0, -t36, 0, 0, 0; 0, 0, -t44, 0, 0, 0; t36, 0, t35, 0, 0, 0; t34, 0, t37, 0, 0, 0; 0, 0, -t43, 0, 0, 0; -qJD(1) * t38, 0, 0, 0, 0, 0; qJD(1) * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t57 = qJ(1) + pkin(10);
	t53 = sin(t57);
	t62 = qJD(1) * t53;
	t55 = cos(t57);
	t61 = qJD(1) * t55;
	t56 = qJ(3) + pkin(11);
	t52 = sin(t56);
	t60 = qJD(3) * t52;
	t54 = cos(t56);
	t59 = qJD(3) * t54;
	t58 = qJD(3) * t55;
	t51 = t53 * t60 - t54 * t61;
	t50 = t52 * t61 + t53 * t59;
	t49 = t52 * t58 + t54 * t62;
	t48 = t52 * t62 - t54 * t58;
	t1 = [t51, 0, t48, 0, 0, 0; -t49, 0, -t50, 0, 0, 0; 0, 0, -t60, 0, 0, 0; t50, 0, t49, 0, 0, 0; t48, 0, t51, 0, 0, 0; 0, 0, -t59, 0, 0, 0; -t62, 0, 0, 0, 0, 0; t61, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (115->15), mult. (54->12), div. (0->0), fcn. (54->4), ass. (0->16)
	t84 = qJ(3) + pkin(11) + qJ(5);
	t80 = sin(t84);
	t85 = qJD(3) + qJD(5);
	t90 = t85 * t80;
	t81 = cos(t84);
	t89 = t85 * t81;
	t86 = qJ(1) + pkin(10);
	t82 = sin(t86);
	t88 = qJD(1) * t82;
	t83 = cos(t86);
	t87 = qJD(1) * t83;
	t79 = -t81 * t87 + t82 * t90;
	t78 = t80 * t87 + t82 * t89;
	t77 = t81 * t88 + t83 * t90;
	t76 = t80 * t88 - t83 * t89;
	t1 = [t79, 0, t76, 0, t76, 0; -t77, 0, -t78, 0, -t78, 0; 0, 0, -t90, 0, -t90, 0; t78, 0, t77, 0, t77, 0; t76, 0, t79, 0, t79, 0; 0, 0, -t89, 0, -t89, 0; -t88, 0, 0, 0, 0, 0; t87, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:16
	% EndTime: 2019-10-10 00:46:16
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (314->29), mult. (226->56), div. (0->0), fcn. (226->6), ass. (0->42)
	t352 = cos(qJ(6));
	t348 = qJ(3) + pkin(11) + qJ(5);
	t345 = cos(t348);
	t355 = qJD(1) * t345 - qJD(6);
	t371 = t352 * t355;
	t356 = qJD(6) * t345 - qJD(1);
	t344 = sin(t348);
	t349 = qJD(3) + qJD(5);
	t351 = sin(qJ(6));
	t368 = t349 * t351;
	t360 = t344 * t368;
	t370 = t356 * t352 - t360;
	t369 = t344 * t349;
	t343 = t349 * t345;
	t367 = t349 * t352;
	t350 = qJ(1) + pkin(10);
	t346 = sin(t350);
	t366 = qJD(1) * t346;
	t347 = cos(t350);
	t365 = qJD(1) * t347;
	t364 = qJD(1) * t351;
	t363 = qJD(1) * t352;
	t362 = qJD(6) * t351;
	t361 = qJD(6) * t352;
	t359 = t344 * t367;
	t358 = t345 * t368;
	t357 = t345 * t367;
	t354 = t355 * t351;
	t353 = t356 * t351 + t359;
	t342 = -t345 * t362 - t359;
	t341 = -t345 * t361 + t360;
	t340 = t345 * t365 - t346 * t369;
	t339 = -t345 * t366 - t347 * t369;
	t338 = -t346 * t357 + (t346 * t362 - t347 * t363) * t344;
	t337 = t346 * t358 + (t346 * t361 + t347 * t364) * t344;
	t336 = -t347 * t357 + (t346 * t363 + t347 * t362) * t344;
	t335 = t347 * t358 + (-t346 * t364 + t347 * t361) * t344;
	t334 = t353 * t346 - t347 * t371;
	t333 = t370 * t346 + t347 * t354;
	t332 = t346 * t371 + t353 * t347;
	t331 = t346 * t354 - t370 * t347;
	t1 = [t334, 0, t336, 0, t336, t331; -t332, 0, t338, 0, t338, -t333; 0, 0, t342, 0, t342, -t344 * t361 - t358; t333, 0, t335, 0, t335, t332; t331, 0, t337, 0, t337, t334; 0, 0, t341, 0, t341, t344 * t362 - t357; -t346 * t343 - t344 * t365, 0, t339, 0, t339, 0; t347 * t343 - t344 * t366, 0, t340, 0, t340, 0; 0, 0, t343, 0, t343, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end