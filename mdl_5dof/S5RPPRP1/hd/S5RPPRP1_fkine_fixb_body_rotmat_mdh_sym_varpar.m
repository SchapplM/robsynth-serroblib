% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:12
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:05
	% EndTime: 2020-11-04 20:12:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:05
	% EndTime: 2020-11-04 20:12:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t43 = cos(qJ(1));
	t42 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t42, t43, 0, 0; -t43, t42, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:05
	% EndTime: 2020-11-04 20:12:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t46 = qJ(1) + pkin(7);
	t45 = cos(t46);
	t44 = sin(t46);
	t1 = [0, 0, 1, qJ(2) + pkin(5) + 0; t44, t45, 0, sin(qJ(1)) * pkin(1) + 0; -t45, t44, 0, -cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:05
	% EndTime: 2020-11-04 20:12:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->13), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t51 = cos(pkin(8));
	t50 = sin(pkin(8));
	t49 = qJ(1) + pkin(7);
	t48 = cos(t49);
	t47 = sin(t49);
	t1 = [t50, t51, 0, qJ(2) + pkin(5) + 0; t47 * t51, -t47 * t50, -t48, t47 * pkin(2) - t48 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; -t48 * t51, t48 * t50, -t47, -t48 * pkin(2) - t47 * qJ(3) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:05
	% EndTime: 2020-11-04 20:12:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (39->20), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t56 = cos(pkin(8));
	t57 = sin(qJ(4));
	t61 = t56 * t57;
	t58 = cos(qJ(4));
	t60 = t56 * t58;
	t55 = sin(pkin(8));
	t59 = pkin(3) * t56 + pkin(6) * t55 + pkin(2);
	t54 = qJ(1) + pkin(7);
	t53 = cos(t54);
	t52 = sin(t54);
	t1 = [t55 * t58, -t55 * t57, -t56, t55 * pkin(3) - t56 * pkin(6) + pkin(5) + qJ(2) + 0; t52 * t60 - t53 * t57, -t52 * t61 - t53 * t58, t52 * t55, sin(qJ(1)) * pkin(1) - t53 * qJ(3) + 0 + t59 * t52; -t52 * t57 - t53 * t60, -t52 * t58 + t53 * t61, -t53 * t55, -cos(qJ(1)) * pkin(1) - t52 * qJ(3) + 0 - t59 * t53; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:05
	% EndTime: 2020-11-04 20:12:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (49->23), mult. (39->26), div. (0->0), fcn. (52->8), ass. (0->14)
	t67 = cos(pkin(8));
	t69 = sin(qJ(4));
	t74 = t67 * t69;
	t70 = cos(qJ(4));
	t73 = t67 * t70;
	t72 = -pkin(4) * t69 - qJ(3);
	t62 = t70 * pkin(4) + pkin(3);
	t66 = sin(pkin(8));
	t68 = -qJ(5) - pkin(6);
	t71 = t62 * t67 - t66 * t68 + pkin(2);
	t65 = qJ(1) + pkin(7);
	t64 = cos(t65);
	t63 = sin(t65);
	t1 = [t66 * t70, -t66 * t69, -t67, t66 * t62 + t67 * t68 + pkin(5) + qJ(2) + 0; t63 * t73 - t64 * t69, -t63 * t74 - t64 * t70, t63 * t66, sin(qJ(1)) * pkin(1) + 0 + t72 * t64 + t71 * t63; -t63 * t69 - t64 * t73, -t63 * t70 + t64 * t74, -t64 * t66, -cos(qJ(1)) * pkin(1) + 0 + t72 * t63 - t71 * t64; 0, 0, 0, 1;];
	Tc_mdh = t1;
end