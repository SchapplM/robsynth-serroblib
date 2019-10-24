% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:42
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRPR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:57
	% EndTime: 2019-10-24 10:42:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:57
	% EndTime: 2019-10-24 10:42:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t4, 0, 0, 0, 0; -t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; t3, 0, 0, 0, 0; -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:57
	% EndTime: 2019-10-24 10:42:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t15 = cos(qJ(1));
	t14 = sin(qJ(1));
	t13 = cos(pkin(8));
	t12 = sin(pkin(8));
	t1 = [0, 0, 0, 0, 0; -t15 * t13, 0, 0, 0, 0; -t14 * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; t15 * t12, 0, 0, 0, 0; t14 * t12, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0; t15, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:57
	% EndTime: 2019-10-24 10:42:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (12->10), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t41 = sin(qJ(3));
	t42 = sin(qJ(1));
	t48 = t42 * t41;
	t43 = cos(qJ(3));
	t47 = t42 * t43;
	t44 = cos(qJ(1));
	t46 = t44 * t41;
	t45 = t44 * t43;
	t40 = cos(pkin(8));
	t39 = sin(pkin(8));
	t38 = -t40 * t45 - t48;
	t37 = t40 * t46 - t47;
	t36 = t40 * t47 - t46;
	t35 = t40 * t48 + t45;
	t1 = [0, 0, -t39 * t41, 0, 0; t38, 0, t35, 0, 0; -t36, 0, -t37, 0, 0; 0, 0, -t39 * t43, 0, 0; t37, 0, t36, 0, 0; t35, 0, t38, 0, 0; 0, 0, 0, 0, 0; -t44 * t39, 0, 0, 0, 0; -t42 * t39, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:57
	% EndTime: 2019-10-24 10:42:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->11), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->16)
	t47 = qJ(3) + pkin(9);
	t45 = sin(t47);
	t50 = sin(qJ(1));
	t55 = t50 * t45;
	t46 = cos(t47);
	t54 = t50 * t46;
	t51 = cos(qJ(1));
	t53 = t51 * t45;
	t52 = t51 * t46;
	t49 = cos(pkin(8));
	t48 = sin(pkin(8));
	t44 = -t49 * t52 - t55;
	t43 = t49 * t53 - t54;
	t42 = t49 * t54 - t53;
	t41 = t49 * t55 + t52;
	t1 = [0, 0, -t48 * t45, 0, 0; t44, 0, t41, 0, 0; -t42, 0, -t43, 0, 0; 0, 0, -t48 * t46, 0, 0; t43, 0, t42, 0, 0; t41, 0, t44, 0, 0; 0, 0, 0, 0, 0; -t51 * t48, 0, 0, 0, 0; -t50 * t48, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:57
	% EndTime: 2019-10-24 10:42:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (74->14), mult. (42->12), div. (0->0), fcn. (72->6), ass. (0->18)
	t63 = qJ(3) + pkin(9) + qJ(5);
	t61 = sin(t63);
	t64 = sin(pkin(8));
	t73 = t64 * t61;
	t62 = cos(t63);
	t72 = t64 * t62;
	t66 = sin(qJ(1));
	t71 = t66 * t61;
	t70 = t66 * t62;
	t67 = cos(qJ(1));
	t69 = t67 * t61;
	t68 = t67 * t62;
	t65 = cos(pkin(8));
	t59 = -t65 * t68 - t71;
	t58 = t65 * t69 - t70;
	t57 = t65 * t70 - t69;
	t56 = t65 * t71 + t68;
	t1 = [0, 0, -t73, 0, -t73; t59, 0, t56, 0, t56; -t57, 0, -t58, 0, -t58; 0, 0, -t72, 0, -t72; t58, 0, t57, 0, t57; t56, 0, t59, 0, t59; 0, 0, 0, 0, 0; -t67 * t64, 0, 0, 0, 0; -t66 * t64, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end