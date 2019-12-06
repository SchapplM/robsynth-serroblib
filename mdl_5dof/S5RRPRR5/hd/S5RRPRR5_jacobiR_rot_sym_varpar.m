% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:48
	% EndTime: 2019-12-05 18:34:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:48
	% EndTime: 2019-12-05 18:34:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t4, 0, 0, 0, 0; -t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; t3, 0, 0, 0, 0; -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:48
	% EndTime: 2019-12-05 18:34:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t9 = qJ(1) + qJ(2);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; -t8, -t8, 0, 0, 0; -t7, -t7, 0, 0, 0; 0, 0, 0, 0, 0; t7, t7, 0, 0, 0; -t8, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:48
	% EndTime: 2019-12-05 18:34:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t27 = qJ(1) + qJ(2);
	t25 = sin(t27);
	t29 = cos(pkin(9));
	t31 = t25 * t29;
	t26 = cos(t27);
	t30 = t26 * t29;
	t28 = sin(pkin(9));
	t24 = t26 * t28;
	t23 = t25 * t28;
	t1 = [0, 0, 0, 0, 0; -t30, -t30, 0, 0, 0; -t31, -t31, 0, 0, 0; 0, 0, 0, 0, 0; t24, t24, 0, 0, 0; t23, t23, 0, 0, 0; 0, 0, 0, 0, 0; -t25, -t25, 0, 0, 0; t26, t26, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:48
	% EndTime: 2019-12-05 18:34:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (39->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t42 = pkin(9) + qJ(4);
	t39 = cos(t42);
	t43 = qJ(1) + qJ(2);
	t40 = sin(t43);
	t45 = t40 * t39;
	t41 = cos(t43);
	t44 = t41 * t39;
	t38 = sin(t42);
	t37 = t41 * t38;
	t36 = t40 * t38;
	t1 = [0, 0, 0, t39, 0; -t44, -t44, 0, t36, 0; -t45, -t45, 0, -t37, 0; 0, 0, 0, -t38, 0; t37, t37, 0, t45, 0; t36, t36, 0, -t44, 0; 0, 0, 0, 0, 0; -t40, -t40, 0, 0, 0; t41, t41, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:48
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (72->14), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->11)
	t52 = pkin(9) + qJ(4) + qJ(5);
	t51 = cos(t52);
	t55 = qJ(1) + qJ(2);
	t54 = cos(t55);
	t56 = t54 * t51;
	t53 = sin(t55);
	t50 = sin(t52);
	t49 = t54 * t50;
	t48 = t53 * t51;
	t47 = t53 * t50;
	t1 = [0, 0, 0, t51, t51; -t56, -t56, 0, t47, t47; -t48, -t48, 0, -t49, -t49; 0, 0, 0, -t50, -t50; t49, t49, 0, t48, t48; t47, t47, 0, -t56, -t56; 0, 0, 0, 0, 0; -t53, -t53, 0, 0, 0; t54, t54, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end