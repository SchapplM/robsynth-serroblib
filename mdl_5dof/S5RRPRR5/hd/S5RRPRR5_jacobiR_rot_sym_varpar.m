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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t4, 0, 0, 0, 0; t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t3, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->3), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t11 = qJ(1) + qJ(2);
	t10 = cos(t11);
	t9 = sin(t11);
	t1 = [0, 0, 0, 0, 0; t10, t10, 0, 0, 0; t9, t9, 0, 0, 0; 0, 0, 0, 0, 0; -t9, -t9, 0, 0, 0; t10, t10, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t29 = qJ(1) + qJ(2);
	t27 = sin(t29);
	t30 = sin(pkin(9));
	t33 = t27 * t30;
	t28 = cos(t29);
	t32 = t28 * t30;
	t31 = cos(pkin(9));
	t26 = t28 * t31;
	t25 = t27 * t31;
	t1 = [0, 0, 0, 0, 0; t26, t26, 0, 0, 0; t25, t25, 0, 0, 0; 0, 0, 0, 0, 0; -t32, -t32, 0, 0, 0; -t33, -t33, 0, 0, 0; 0, 0, 0, 0, 0; t27, t27, 0, 0, 0; -t28, -t28, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (39->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t43 = pkin(9) + qJ(4);
	t39 = sin(t43);
	t44 = qJ(1) + qJ(2);
	t41 = sin(t44);
	t46 = t41 * t39;
	t42 = cos(t44);
	t45 = t42 * t39;
	t40 = cos(t43);
	t38 = t42 * t40;
	t37 = t41 * t40;
	t1 = [0, 0, 0, t40, 0; t38, t38, 0, -t46, 0; t37, t37, 0, t45, 0; 0, 0, 0, -t39, 0; -t45, -t45, 0, -t37, 0; -t46, -t46, 0, t38, 0; 0, 0, 0, 0, 0; t41, t41, 0, 0, 0; -t42, -t42, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (72->14), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->11)
	t53 = pkin(9) + qJ(4) + qJ(5);
	t51 = sin(t53);
	t56 = qJ(1) + qJ(2);
	t54 = sin(t56);
	t57 = t54 * t51;
	t55 = cos(t56);
	t52 = cos(t53);
	t50 = t55 * t52;
	t49 = t55 * t51;
	t48 = t54 * t52;
	t1 = [0, 0, 0, t52, t52; t50, t50, 0, -t57, -t57; t48, t48, 0, t49, t49; 0, 0, 0, -t51, -t51; -t49, -t49, 0, -t48, -t48; -t57, -t57, 0, t50, t50; 0, 0, 0, 0, 0; t54, t54, 0, 0, 0; -t55, -t55, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end