% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPPRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:52
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:59
	% EndTime: 2020-11-04 19:51:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:59
	% EndTime: 2020-11-04 19:51:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t52 = cos(pkin(7));
	t51 = sin(pkin(7));
	t1 = [t52, -t51, 0, 0; t51, t52, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:59
	% EndTime: 2020-11-04 19:51:59
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t56 = cos(pkin(7));
	t55 = cos(pkin(8));
	t54 = sin(pkin(7));
	t53 = sin(pkin(8));
	t1 = [t56 * t55, -t56 * t53, t54, pkin(1) * t56 + qJ(2) * t54 + 0; t54 * t55, -t54 * t53, -t56, pkin(1) * t54 - qJ(2) * t56 + 0; t53, t55, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:59
	% EndTime: 2020-11-04 19:51:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t60 = sin(pkin(7));
	t62 = cos(pkin(8));
	t66 = t60 * t62;
	t58 = sin(pkin(9));
	t63 = cos(pkin(7));
	t65 = t63 * t58;
	t61 = cos(pkin(9));
	t64 = t63 * t61;
	t59 = sin(pkin(8));
	t57 = t62 * pkin(2) + t59 * qJ(3) + pkin(1);
	t1 = [t60 * t58 + t62 * t64, t60 * t61 - t62 * t65, t63 * t59, t60 * qJ(2) + t57 * t63 + 0; t61 * t66 - t65, -t58 * t66 - t64, t60 * t59, -t63 * qJ(2) + t57 * t60 + 0; t59 * t61, -t59 * t58, -t62, t59 * pkin(2) - t62 * qJ(3) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:59
	% EndTime: 2020-11-04 19:51:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t74 = sin(pkin(7));
	t75 = cos(pkin(8));
	t80 = t74 * t75;
	t72 = pkin(9) + qJ(4);
	t70 = sin(t72);
	t76 = cos(pkin(7));
	t79 = t76 * t70;
	t71 = cos(t72);
	t78 = t76 * t71;
	t77 = qJ(3) + pkin(5);
	t73 = sin(pkin(8));
	t69 = cos(pkin(9)) * pkin(3) + pkin(2);
	t68 = sin(pkin(9)) * pkin(3) + qJ(2);
	t67 = t69 * t75 + t77 * t73 + pkin(1);
	t1 = [t74 * t70 + t75 * t78, t74 * t71 - t75 * t79, t76 * t73, t67 * t76 + t68 * t74 + 0; t71 * t80 - t79, -t70 * t80 - t78, t74 * t73, t67 * t74 - t68 * t76 + 0; t73 * t71, -t73 * t70, -t75, t73 * t69 - t75 * t77 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:59
	% EndTime: 2020-11-04 19:51:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (66->29), mult. (72->40), div. (0->0), fcn. (97->10), ass. (0->23)
	t91 = sin(pkin(8));
	t96 = sin(qJ(5));
	t102 = t91 * t96;
	t97 = cos(qJ(5));
	t101 = t91 * t97;
	t92 = sin(pkin(7));
	t93 = cos(pkin(8));
	t100 = t92 * t93;
	t90 = pkin(9) + qJ(4);
	t88 = sin(t90);
	t94 = cos(pkin(7));
	t99 = t94 * t88;
	t89 = cos(t90);
	t98 = t94 * t89;
	t95 = qJ(3) + pkin(5);
	t87 = cos(pkin(9)) * pkin(3) + pkin(2);
	t86 = sin(pkin(9)) * pkin(3) + qJ(2);
	t85 = t87 * t93 + t95 * t91 + pkin(1);
	t84 = t92 * t88 + t93 * t98;
	t83 = -t92 * t89 + t93 * t99;
	t82 = t89 * t100 - t99;
	t81 = t88 * t100 + t98;
	t1 = [t94 * t102 + t84 * t97, t94 * t101 - t84 * t96, t83, t84 * pkin(4) + t83 * pkin(6) + t85 * t94 + t86 * t92 + 0; t92 * t102 + t82 * t97, t92 * t101 - t82 * t96, t81, t82 * pkin(4) + t81 * pkin(6) + t85 * t92 - t86 * t94 + 0; t89 * t101 - t93 * t96, -t89 * t102 - t93 * t97, t91 * t88, -t93 * t95 + qJ(1) + 0 + (pkin(4) * t89 + pkin(6) * t88 + t87) * t91; 0, 0, 0, 1;];
	Tc_mdh = t1;
end