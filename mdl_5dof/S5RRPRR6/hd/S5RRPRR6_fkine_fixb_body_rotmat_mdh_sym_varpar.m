% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:16:42
	% EndTime: 2022-01-20 11:16:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:16:42
	% EndTime: 2022-01-20 11:16:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t1 = [t47, -t46, 0, 0; t46, t47, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:16:42
	% EndTime: 2022-01-20 11:16:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t50 = qJ(1) + qJ(2);
	t49 = cos(t50);
	t48 = sin(t50);
	t1 = [t49, -t48, 0, cos(qJ(1)) * pkin(1) + 0; t48, t49, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:16:42
	% EndTime: 2022-01-20 11:16:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t55 = cos(pkin(9));
	t54 = sin(pkin(9));
	t53 = qJ(1) + qJ(2);
	t52 = cos(t53);
	t51 = sin(t53);
	t1 = [t52 * t55, -t52 * t54, t51, t52 * pkin(2) + t51 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; t51 * t55, -t51 * t54, -t52, t51 * pkin(2) - t52 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; t54, t55, 0, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:16:42
	% EndTime: 2022-01-20 11:16:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t60 = cos(pkin(9));
	t61 = sin(qJ(4));
	t65 = t60 * t61;
	t62 = cos(qJ(4));
	t64 = t60 * t62;
	t59 = sin(pkin(9));
	t63 = pkin(3) * t60 + pkin(7) * t59 + pkin(2);
	t58 = qJ(1) + qJ(2);
	t57 = cos(t58);
	t56 = sin(t58);
	t1 = [t56 * t61 + t57 * t64, t56 * t62 - t57 * t65, t57 * t59, cos(qJ(1)) * pkin(1) + t56 * qJ(3) + 0 + t63 * t57; t56 * t64 - t57 * t61, -t56 * t65 - t57 * t62, t56 * t59, sin(qJ(1)) * pkin(1) - t57 * qJ(3) + 0 + t63 * t56; t59 * t62, -t59 * t61, -t60, t59 * pkin(3) - t60 * pkin(7) + pkin(5) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:16:42
	% EndTime: 2022-01-20 11:16:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->23), mult. (39->26), div. (0->0), fcn. (52->10), ass. (0->15)
	t72 = qJ(1) + qJ(2);
	t68 = sin(t72);
	t74 = cos(pkin(9));
	t80 = t68 * t74;
	t70 = cos(t72);
	t79 = t70 * t74;
	t78 = pkin(4) * sin(qJ(4)) + qJ(3);
	t66 = cos(qJ(4)) * pkin(4) + pkin(3);
	t73 = sin(pkin(9));
	t76 = -pkin(8) - pkin(7);
	t77 = t66 * t74 - t73 * t76 + pkin(2);
	t71 = qJ(4) + qJ(5);
	t69 = cos(t71);
	t67 = sin(t71);
	t1 = [t68 * t67 + t69 * t79, -t67 * t79 + t68 * t69, t70 * t73, cos(qJ(1)) * pkin(1) + 0 + t78 * t68 + t77 * t70; -t70 * t67 + t69 * t80, -t67 * t80 - t70 * t69, t68 * t73, sin(qJ(1)) * pkin(1) + 0 - t78 * t70 + t77 * t68; t73 * t69, -t73 * t67, -t74, t73 * t66 + t74 * t76 + pkin(5) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end