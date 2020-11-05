% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP13 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:25
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRP13_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP13_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:52
	% EndTime: 2020-11-04 20:25:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:53
	% EndTime: 2020-11-04 20:25:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t43 = cos(qJ(1));
	t42 = sin(qJ(1));
	t1 = [t43, -t42, 0, 0; t42, t43, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:53
	% EndTime: 2020-11-04 20:25:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t45 = cos(qJ(1));
	t44 = sin(qJ(1));
	t1 = [0, -t45, t44, t45 * pkin(1) + t44 * qJ(2) + 0; 0, -t44, -t45, t44 * pkin(1) - t45 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:53
	% EndTime: 2020-11-04 20:25:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t50 = pkin(1) + pkin(6);
	t49 = cos(qJ(1));
	t48 = cos(qJ(3));
	t47 = sin(qJ(1));
	t46 = sin(qJ(3));
	t1 = [t47 * t46, t47 * t48, t49, t47 * qJ(2) + t50 * t49 + 0; -t49 * t46, -t49 * t48, t47, -t49 * qJ(2) + t50 * t47 + 0; t48, -t46, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:53
	% EndTime: 2020-11-04 20:25:53
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t52 = sin(qJ(4));
	t54 = sin(qJ(1));
	t62 = t54 * t52;
	t55 = cos(qJ(4));
	t61 = t54 * t55;
	t57 = cos(qJ(1));
	t60 = t57 * t52;
	t59 = t57 * t55;
	t58 = pkin(1) + pkin(6);
	t56 = cos(qJ(3));
	t53 = sin(qJ(3));
	t51 = -t53 * pkin(3) + t56 * pkin(7) - qJ(2);
	t1 = [t53 * t61 + t60, -t53 * t62 + t59, -t54 * t56, -t51 * t54 + t58 * t57 + 0; -t53 * t59 + t62, t53 * t60 + t61, t57 * t56, t51 * t57 + t58 * t54 + 0; t56 * t55, -t56 * t52, t53, t56 * pkin(3) + t53 * pkin(7) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:53
	% EndTime: 2020-11-04 20:25:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (29->20), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->14)
	t67 = sin(qJ(4));
	t69 = sin(qJ(1));
	t77 = t69 * t67;
	t70 = cos(qJ(4));
	t76 = t69 * t70;
	t72 = cos(qJ(1));
	t75 = t72 * t67;
	t74 = t72 * t70;
	t64 = -t70 * pkin(4) - t67 * qJ(5) - pkin(3);
	t68 = sin(qJ(3));
	t71 = cos(qJ(3));
	t73 = t71 * pkin(7) + t64 * t68 - qJ(2);
	t63 = t67 * pkin(4) - qJ(5) * t70 + pkin(1) + pkin(6);
	t1 = [t68 * t76 + t75, -t69 * t71, t68 * t77 - t74, t63 * t72 - t73 * t69 + 0; -t68 * t74 + t77, t72 * t71, -t68 * t75 - t76, t63 * t69 + t73 * t72 + 0; t71 * t70, t68, t71 * t67, t68 * pkin(7) - t64 * t71 + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end