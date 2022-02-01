% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:16:08
	% EndTime: 2022-01-23 09:16:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:16:08
	% EndTime: 2022-01-23 09:16:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [t46, -t45, 0, 0; t45, t46, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:16:08
	% EndTime: 2022-01-23 09:16:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t50 = cos(qJ(1));
	t49 = sin(qJ(1));
	t48 = cos(pkin(8));
	t47 = sin(pkin(8));
	t1 = [t50 * t48, -t50 * t47, t49, pkin(1) * t50 + qJ(2) * t49 + 0; t49 * t48, -t49 * t47, -t50, pkin(1) * t49 - qJ(2) * t50 + 0; t47, t48, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:16:08
	% EndTime: 2022-01-23 09:16:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->15), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->12)
	t52 = sin(pkin(9));
	t56 = sin(qJ(1));
	t61 = t56 * t52;
	t54 = cos(pkin(9));
	t60 = t56 * t54;
	t57 = cos(qJ(1));
	t59 = t57 * t52;
	t58 = t57 * t54;
	t55 = cos(pkin(8));
	t53 = sin(pkin(8));
	t51 = pkin(2) * t55 + t53 * qJ(3) + pkin(1);
	t1 = [t55 * t58 + t61, -t55 * t59 + t60, t57 * t53, t56 * qJ(2) + t51 * t57 + 0; t55 * t60 - t59, -t55 * t61 - t58, t56 * t53, -t57 * qJ(2) + t51 * t56 + 0; t53 * t54, -t53 * t52, -t55, t53 * pkin(2) - t55 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:16:08
	% EndTime: 2022-01-23 09:16:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t67 = pkin(9) + qJ(4);
	t65 = sin(t67);
	t71 = sin(qJ(1));
	t76 = t71 * t65;
	t66 = cos(t67);
	t75 = t71 * t66;
	t72 = cos(qJ(1));
	t74 = t72 * t65;
	t73 = t72 * t66;
	t70 = qJ(3) + pkin(6);
	t69 = cos(pkin(8));
	t68 = sin(pkin(8));
	t64 = cos(pkin(9)) * pkin(3) + pkin(2);
	t63 = sin(pkin(9)) * pkin(3) + qJ(2);
	t62 = t64 * t69 + t70 * t68 + pkin(1);
	t1 = [t69 * t73 + t76, -t69 * t74 + t75, t72 * t68, t62 * t72 + t63 * t71 + 0; t69 * t75 - t74, -t69 * t76 - t73, t71 * t68, t62 * t71 - t63 * t72 + 0; t68 * t66, -t68 * t65, -t69, t68 * t64 - t69 * t70 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:16:08
	% EndTime: 2022-01-23 09:16:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->22), mult. (42->24), div. (0->0), fcn. (55->10), ass. (0->17)
	t83 = pkin(9) + qJ(4);
	t81 = qJ(5) + t83;
	t79 = sin(t81);
	t86 = sin(qJ(1));
	t93 = t86 * t79;
	t80 = cos(t81);
	t92 = t86 * t80;
	t87 = cos(qJ(1));
	t91 = t87 * t79;
	t90 = t87 * t80;
	t89 = qJ(2) + pkin(4) * sin(t83) + sin(pkin(9)) * pkin(3);
	t77 = pkin(4) * cos(t83) + cos(pkin(9)) * pkin(3) + pkin(2);
	t82 = -pkin(7) - pkin(6) - qJ(3);
	t84 = sin(pkin(8));
	t85 = cos(pkin(8));
	t88 = t77 * t85 - t82 * t84 + pkin(1);
	t1 = [t85 * t90 + t93, -t85 * t91 + t92, t87 * t84, t89 * t86 + t88 * t87 + 0; t85 * t92 - t91, -t85 * t93 - t90, t86 * t84, t88 * t86 - t89 * t87 + 0; t84 * t80, -t84 * t79, -t85, t84 * t77 + t85 * t82 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end