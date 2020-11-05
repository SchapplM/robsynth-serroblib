% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:11
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:32
	% EndTime: 2020-11-04 20:11:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:32
	% EndTime: 2020-11-04 20:11:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(qJ(1));
	t35 = sin(qJ(1));
	t1 = [t36, -t35, 0, 0; t35, t36, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:32
	% EndTime: 2020-11-04 20:11:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [t38, 0, t37, t38 * pkin(1) + t37 * qJ(2) + 0; t37, 0, -t38, t37 * pkin(1) - t38 * qJ(2) + 0; 0, 1, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:32
	% EndTime: 2020-11-04 20:11:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t45 = pkin(1) + pkin(2);
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t42 = cos(pkin(7));
	t41 = sin(pkin(7));
	t40 = -t44 * t41 + t43 * t42;
	t39 = -t43 * t41 - t44 * t42;
	t1 = [-t39, t40, 0, t43 * qJ(2) + t45 * t44 + 0; t40, t39, 0, -t44 * qJ(2) + t45 * t43 + 0; 0, 0, -1, -qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:32
	% EndTime: 2020-11-04 20:11:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (26->17), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->11)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t54 = cos(pkin(7));
	t53 = cos(pkin(8));
	t52 = sin(pkin(7));
	t51 = sin(pkin(8));
	t49 = -pkin(3) * t52 + qJ(4) * t54 - qJ(2);
	t48 = pkin(3) * t54 + qJ(4) * t52 + pkin(1) + pkin(2);
	t47 = t55 * t52 + t54 * t56;
	t46 = t52 * t56 - t55 * t54;
	t1 = [t47 * t53, -t47 * t51, t46, t48 * t56 - t49 * t55 + 0; -t46 * t53, t46 * t51, t47, t48 * t55 + t49 * t56 + 0; -t51, -t53, 0, -qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:32
	% EndTime: 2020-11-04 20:11:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (41->21), mult. (33->18), div. (0->0), fcn. (47->8), ass. (0->14)
	t61 = cos(pkin(8)) * pkin(4) + pkin(3);
	t65 = sin(pkin(7));
	t66 = cos(pkin(7));
	t67 = qJ(4) + pkin(6);
	t70 = t61 * t65 - t67 * t66 + qJ(2);
	t69 = cos(qJ(1));
	t68 = sin(qJ(1));
	t64 = pkin(8) + qJ(5);
	t63 = cos(t64);
	t62 = sin(t64);
	t59 = t68 * t65 + t69 * t66;
	t58 = t69 * t65 - t68 * t66;
	t57 = t61 * t66 + t67 * t65 + pkin(1) + pkin(2);
	t1 = [t59 * t63, -t59 * t62, t58, t57 * t69 + t70 * t68 + 0; -t58 * t63, t58 * t62, t59, t57 * t68 - t70 * t69 + 0; -t62, -t63, 0, -sin(pkin(8)) * pkin(4) - qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end