% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:00
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:33
	% EndTime: 2020-11-04 20:00:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:33
	% EndTime: 2020-11-04 20:00:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t47 = cos(pkin(8));
	t46 = sin(pkin(8));
	t1 = [t47, -t46, 0, 0; t46, t47, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:33
	% EndTime: 2020-11-04 20:00:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t51 = cos(qJ(2));
	t50 = sin(qJ(2));
	t49 = cos(pkin(8));
	t48 = sin(pkin(8));
	t1 = [t49 * t51, -t49 * t50, t48, t49 * pkin(1) + t48 * pkin(5) + 0; t48 * t51, -t48 * t50, -t49, t48 * pkin(1) - t49 * pkin(5) + 0; t50, t51, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:33
	% EndTime: 2020-11-04 20:00:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->10)
	t53 = sin(pkin(8));
	t57 = cos(qJ(2));
	t60 = t53 * t57;
	t55 = cos(pkin(8));
	t59 = t55 * t57;
	t56 = sin(qJ(2));
	t58 = pkin(2) * t57 + qJ(3) * t56 + pkin(1);
	t54 = cos(pkin(9));
	t52 = sin(pkin(9));
	t1 = [t53 * t52 + t54 * t59, -t52 * t59 + t53 * t54, t55 * t56, t53 * pkin(5) + t58 * t55 + 0; -t55 * t52 + t54 * t60, -t52 * t60 - t55 * t54, t53 * t56, -t55 * pkin(5) + t58 * t53 + 0; t56 * t54, -t56 * t52, -t57, t56 * pkin(2) - t57 * qJ(3) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:33
	% EndTime: 2020-11-04 20:00:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (35->24), div. (0->0), fcn. (48->8), ass. (0->14)
	t66 = sin(pkin(8));
	t70 = cos(qJ(2));
	t73 = t66 * t70;
	t67 = cos(pkin(8));
	t72 = t67 * t70;
	t62 = cos(pkin(9)) * pkin(3) + pkin(2);
	t68 = qJ(3) + pkin(6);
	t69 = sin(qJ(2));
	t71 = t62 * t70 + t68 * t69 + pkin(1);
	t65 = pkin(9) + qJ(4);
	t64 = cos(t65);
	t63 = sin(t65);
	t61 = sin(pkin(9)) * pkin(3) + pkin(5);
	t1 = [t66 * t63 + t64 * t72, -t63 * t72 + t66 * t64, t67 * t69, t61 * t66 + t71 * t67 + 0; -t67 * t63 + t64 * t73, -t63 * t73 - t67 * t64, t66 * t69, -t61 * t67 + t71 * t66 + 0; t69 * t64, -t69 * t63, -t70, t69 * t62 - t70 * t68 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:33
	% EndTime: 2020-11-04 20:00:33
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->22), mult. (42->26), div. (0->0), fcn. (55->10), ass. (0->15)
	t80 = pkin(9) + qJ(4);
	t88 = pkin(5) + pkin(4) * sin(t80) + sin(pkin(9)) * pkin(3);
	t81 = sin(pkin(8));
	t84 = cos(qJ(2));
	t87 = t81 * t84;
	t82 = cos(pkin(8));
	t86 = t82 * t84;
	t74 = pkin(4) * cos(t80) + cos(pkin(9)) * pkin(3) + pkin(2);
	t79 = -pkin(7) - pkin(6) - qJ(3);
	t83 = sin(qJ(2));
	t85 = t74 * t84 - t79 * t83 + pkin(1);
	t78 = qJ(5) + t80;
	t77 = cos(t78);
	t76 = sin(t78);
	t1 = [t81 * t76 + t77 * t86, -t76 * t86 + t81 * t77, t82 * t83, t88 * t81 + t85 * t82 + 0; -t82 * t76 + t77 * t87, -t76 * t87 - t82 * t77, t81 * t83, t85 * t81 - t88 * t82 + 0; t83 * t77, -t83 * t76, -t84, t83 * t74 + t84 * t79 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end