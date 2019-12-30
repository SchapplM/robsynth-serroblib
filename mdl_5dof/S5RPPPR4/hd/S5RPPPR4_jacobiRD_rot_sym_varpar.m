% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:47
	% EndTime: 2019-12-29 15:48:47
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:47
	% EndTime: 2019-12-29 15:48:47
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:47
	% EndTime: 2019-12-29 15:48:47
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(7);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0; -t37, 0, 0, 0, 0; 0, 0, 0, 0, 0; t37, 0, 0, 0, 0; -t36, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:43
	% EndTime: 2019-12-29 15:48:43
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t16 = qJ(1) + pkin(7);
	t17 = qJD(1) * sin(t16);
	t14 = qJD(1) * cos(t16);
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t14, 0, 0, 0, 0; t17, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t17, 0, 0, 0, 0; t14, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:52
	% EndTime: 2019-12-29 15:48:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->5), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t28 = qJD(1) * sin(pkin(8));
	t27 = qJD(1) * cos(pkin(8));
	t24 = qJ(1) + pkin(7);
	t23 = cos(t24);
	t22 = sin(t24);
	t1 = [-t22 * t28, 0, 0, 0, 0; t23 * t28, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t22 * t27, 0, 0, 0, 0; t23 * t27, 0, 0, 0, 0; 0, 0, 0, 0, 0; -qJD(1) * t23, 0, 0, 0, 0; -qJD(1) * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:42
	% EndTime: 2019-12-29 15:48:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t51 = qJ(1) + pkin(7);
	t47 = sin(t51);
	t56 = qJD(1) * t47;
	t49 = cos(t51);
	t55 = qJD(1) * t49;
	t50 = pkin(8) + qJ(5);
	t46 = sin(t50);
	t54 = qJD(5) * t46;
	t48 = cos(t50);
	t53 = qJD(5) * t48;
	t52 = qJD(5) * t49;
	t45 = -t47 * t54 + t48 * t55;
	t44 = t46 * t55 + t47 * t53;
	t43 = t46 * t52 + t48 * t56;
	t42 = -t46 * t56 + t48 * t52;
	t1 = [t42, 0, 0, 0, t45; t44, 0, 0, 0, t43; 0, 0, 0, 0, -t53; -t43, 0, 0, 0, -t44; t45, 0, 0, 0, t42; 0, 0, 0, 0, t54; -t55, 0, 0, 0, 0; -t56, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end