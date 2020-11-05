% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:10
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:27
	% EndTime: 2020-11-04 20:10:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:27
	% EndTime: 2020-11-04 20:10:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t43, t44, 0, 0; -t44, t43, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:27
	% EndTime: 2020-11-04 20:10:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t47 = qJ(1) + pkin(7);
	t46 = cos(t47);
	t45 = sin(t47);
	t1 = [0, 0, 1, qJ(2) + pkin(5) + 0; t45, t46, 0, sin(qJ(1)) * pkin(1) + 0; -t46, t45, 0, -cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:27
	% EndTime: 2020-11-04 20:10:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->13), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t52 = cos(pkin(8));
	t51 = sin(pkin(8));
	t50 = qJ(1) + pkin(7);
	t49 = cos(t50);
	t48 = sin(t50);
	t1 = [t51, t52, 0, qJ(2) + pkin(5) + 0; t48 * t52, -t48 * t51, -t49, t48 * pkin(2) - t49 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; -t49 * t52, t49 * t51, -t48, -t49 * pkin(2) - t48 * qJ(3) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:27
	% EndTime: 2020-11-04 20:10:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (39->20), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t56 = sin(pkin(9));
	t59 = cos(pkin(8));
	t62 = t56 * t59;
	t58 = cos(pkin(9));
	t61 = t58 * t59;
	t57 = sin(pkin(8));
	t60 = pkin(3) * t59 + qJ(4) * t57 + pkin(2);
	t55 = qJ(1) + pkin(7);
	t54 = cos(t55);
	t53 = sin(t55);
	t1 = [t57 * t58, -t57 * t56, -t59, t57 * pkin(3) - t59 * qJ(4) + pkin(5) + qJ(2) + 0; t53 * t61 - t54 * t56, -t53 * t62 - t54 * t58, t53 * t57, sin(qJ(1)) * pkin(1) - t54 * qJ(3) + 0 + t60 * t53; -t53 * t56 - t54 * t61, -t53 * t58 + t54 * t62, -t54 * t57, -cos(qJ(1)) * pkin(1) - t53 * qJ(3) + 0 - t60 * t54; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:27
	% EndTime: 2020-11-04 20:10:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (59->24), mult. (39->26), div. (0->0), fcn. (52->10), ass. (0->15)
	t69 = qJ(1) + pkin(7);
	t65 = sin(t69);
	t72 = cos(pkin(8));
	t77 = t65 * t72;
	t67 = cos(t69);
	t76 = t67 * t72;
	t75 = -sin(pkin(9)) * pkin(4) - qJ(3);
	t63 = cos(pkin(9)) * pkin(4) + pkin(3);
	t71 = sin(pkin(8));
	t73 = -pkin(6) - qJ(4);
	t74 = t63 * t72 - t71 * t73 + pkin(2);
	t68 = pkin(9) + qJ(5);
	t66 = cos(t68);
	t64 = sin(t68);
	t1 = [t71 * t66, -t71 * t64, -t72, t71 * t63 + t72 * t73 + pkin(5) + qJ(2) + 0; -t67 * t64 + t66 * t77, -t64 * t77 - t67 * t66, t65 * t71, sin(qJ(1)) * pkin(1) + 0 + t75 * t67 + t74 * t65; -t65 * t64 - t66 * t76, t64 * t76 - t65 * t66, -t67 * t71, -cos(qJ(1)) * pkin(1) + 0 + t75 * t65 - t74 * t67; 0, 0, 0, 1;];
	Tc_mdh = t1;
end