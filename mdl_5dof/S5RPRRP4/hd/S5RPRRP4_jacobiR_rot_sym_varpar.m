% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRP4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
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
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t15 = cos(qJ(1));
	t14 = sin(qJ(1));
	t13 = cos(pkin(8));
	t12 = sin(pkin(8));
	t1 = [0, 0, 0, 0, 0; t15 * t13, 0, 0, 0, 0; t14 * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t15 * t12, 0, 0, 0, 0; -t14 * t12, 0, 0, 0, 0; 0, 0, 0, 0, 0; t14, 0, 0, 0, 0; -t15, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->8), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
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
	t38 = t40 * t45 + t48;
	t37 = t40 * t46 - t47;
	t36 = t40 * t47 - t46;
	t35 = -t40 * t48 - t45;
	t1 = [0, 0, -t39 * t41, 0, 0; t38, 0, t35, 0, 0; t36, 0, t37, 0, 0; 0, 0, -t39 * t43, 0, 0; -t37, 0, -t36, 0, 0; t35, 0, t38, 0, 0; 0, 0, 0, 0, 0; t44 * t39, 0, 0, 0, 0; t42 * t39, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->12), mult. (42->12), div. (0->0), fcn. (72->6), ass. (0->18)
	t60 = qJ(3) + qJ(4);
	t58 = sin(t60);
	t61 = sin(pkin(8));
	t70 = t61 * t58;
	t59 = cos(t60);
	t69 = t61 * t59;
	t63 = sin(qJ(1));
	t68 = t63 * t58;
	t67 = t63 * t59;
	t64 = cos(qJ(1));
	t66 = t64 * t58;
	t65 = t64 * t59;
	t62 = cos(pkin(8));
	t56 = t62 * t65 + t68;
	t55 = t62 * t66 - t67;
	t54 = t62 * t67 - t66;
	t53 = -t62 * t68 - t65;
	t1 = [0, 0, -t70, -t70, 0; t56, 0, t53, t53, 0; t54, 0, t55, t55, 0; 0, 0, -t69, -t69, 0; -t55, 0, -t54, -t54, 0; t53, 0, t56, t56, 0; 0, 0, 0, 0, 0; t64 * t61, 0, 0, 0, 0; t63 * t61, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (44->12), mult. (42->12), div. (0->0), fcn. (72->6), ass. (0->18)
	t62 = qJ(3) + qJ(4);
	t60 = sin(t62);
	t63 = sin(pkin(8));
	t72 = t63 * t60;
	t61 = cos(t62);
	t71 = t63 * t61;
	t65 = sin(qJ(1));
	t70 = t65 * t60;
	t69 = t65 * t61;
	t66 = cos(qJ(1));
	t68 = t66 * t60;
	t67 = t66 * t61;
	t64 = cos(pkin(8));
	t58 = t64 * t67 + t70;
	t57 = t64 * t68 - t69;
	t56 = t64 * t69 - t68;
	t55 = -t64 * t70 - t67;
	t1 = [0, 0, -t72, -t72, 0; t58, 0, t55, t55, 0; t56, 0, t57, t57, 0; 0, 0, -t71, -t71, 0; -t57, 0, -t56, -t56, 0; t55, 0, t58, t58, 0; 0, 0, 0, 0, 0; t66 * t63, 0, 0, 0, 0; t65 * t63, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end