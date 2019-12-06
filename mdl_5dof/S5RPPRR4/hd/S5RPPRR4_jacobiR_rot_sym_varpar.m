% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:45:55
	% EndTime: 2019-12-05 17:45:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:45:55
	% EndTime: 2019-12-05 17:45:55
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
	% StartTime: 2019-12-05 17:45:55
	% EndTime: 2019-12-05 17:45:55
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
	% StartTime: 2019-12-05 17:45:55
	% EndTime: 2019-12-05 17:45:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t22 = sin(pkin(9));
	t26 = sin(qJ(1));
	t31 = t26 * t22;
	t24 = cos(pkin(9));
	t30 = t26 * t24;
	t27 = cos(qJ(1));
	t29 = t27 * t22;
	t28 = t27 * t24;
	t25 = cos(pkin(8));
	t23 = sin(pkin(8));
	t1 = [0, 0, 0, 0, 0; -t25 * t28 - t31, 0, 0, 0, 0; -t25 * t30 + t29, 0, 0, 0, 0; 0, 0, 0, 0, 0; t25 * t29 - t30, 0, 0, 0, 0; t25 * t31 + t28, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t27 * t23, 0, 0, 0, 0; -t26 * t23, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:45:55
	% EndTime: 2019-12-05 17:45:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->11), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->16)
	t45 = pkin(9) + qJ(4);
	t43 = sin(t45);
	t48 = sin(qJ(1));
	t53 = t48 * t43;
	t44 = cos(t45);
	t52 = t48 * t44;
	t49 = cos(qJ(1));
	t51 = t49 * t43;
	t50 = t49 * t44;
	t47 = cos(pkin(8));
	t46 = sin(pkin(8));
	t42 = -t47 * t50 - t53;
	t41 = t47 * t51 - t52;
	t40 = t47 * t52 - t51;
	t39 = t47 * t53 + t50;
	t1 = [0, 0, 0, -t46 * t43, 0; t42, 0, 0, t39, 0; -t40, 0, 0, -t41, 0; 0, 0, 0, -t46 * t44, 0; t41, 0, 0, t40, 0; t39, 0, 0, t42, 0; 0, 0, 0, 0, 0; -t49 * t46, 0, 0, 0, 0; -t48 * t46, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:45:55
	% EndTime: 2019-12-05 17:45:55
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (74->14), mult. (42->12), div. (0->0), fcn. (72->6), ass. (0->18)
	t63 = pkin(9) + qJ(4) + qJ(5);
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
	t1 = [0, 0, 0, -t73, -t73; t59, 0, 0, t56, t56; -t57, 0, 0, -t58, -t58; 0, 0, 0, -t72, -t72; t58, 0, 0, t57, t57; t56, 0, 0, t59, t59; 0, 0, 0, 0, 0; -t67 * t64, 0, 0, 0, 0; -t66 * t64, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end