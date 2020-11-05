% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:25
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:20
	% EndTime: 2020-11-04 20:25:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:20
	% EndTime: 2020-11-04 20:25:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t1 = [t51, -t50, 0, 0; t50, t51, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:20
	% EndTime: 2020-11-04 20:25:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t55 = cos(qJ(1));
	t54 = sin(qJ(1));
	t53 = cos(pkin(8));
	t52 = sin(pkin(8));
	t1 = [t55 * t53, -t55 * t52, t54, t55 * pkin(1) + t54 * qJ(2) + 0; t54 * t53, -t54 * t52, -t55, t54 * pkin(1) - t55 * qJ(2) + 0; t52, t53, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:20
	% EndTime: 2020-11-04 20:25:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t60 = pkin(6) + qJ(2);
	t59 = pkin(8) + qJ(3);
	t58 = cos(t59);
	t57 = sin(t59);
	t56 = cos(pkin(8)) * pkin(2) + pkin(1);
	t1 = [t62 * t58, -t62 * t57, t61, t56 * t62 + t60 * t61 + 0; t61 * t58, -t61 * t57, -t62, t56 * t61 - t60 * t62 + 0; t57, t58, 0, sin(pkin(8)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:20
	% EndTime: 2020-11-04 20:25:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t68 = sin(qJ(4));
	t69 = sin(qJ(1));
	t76 = t69 * t68;
	t70 = cos(qJ(4));
	t75 = t69 * t70;
	t71 = cos(qJ(1));
	t74 = t71 * t68;
	t73 = t71 * t70;
	t66 = pkin(8) + qJ(3);
	t64 = sin(t66);
	t65 = cos(t66);
	t72 = pkin(3) * t65 + pkin(7) * t64 + cos(pkin(8)) * pkin(2) + pkin(1);
	t67 = pkin(6) + qJ(2);
	t1 = [t65 * t73 + t76, -t65 * t74 + t75, t71 * t64, t67 * t69 + t72 * t71 + 0; t65 * t75 - t74, -t65 * t76 - t73, t69 * t64, -t71 * t67 + t72 * t69 + 0; t64 * t70, -t64 * t68, -t65, t64 * pkin(3) - t65 * pkin(7) + sin(pkin(8)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:20
	% EndTime: 2020-11-04 20:25:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (52->24), mult. (53->28), div. (0->0), fcn. (70->8), ass. (0->18)
	t86 = sin(qJ(4));
	t87 = sin(qJ(1));
	t94 = t87 * t86;
	t88 = cos(qJ(4));
	t93 = t87 * t88;
	t89 = cos(qJ(1));
	t92 = t89 * t86;
	t91 = t89 * t88;
	t84 = pkin(8) + qJ(3);
	t82 = sin(t84);
	t83 = cos(t84);
	t90 = pkin(3) * t83 + pkin(7) * t82 + cos(pkin(8)) * pkin(2) + pkin(1);
	t85 = pkin(6) + qJ(2);
	t80 = t83 * t91 + t94;
	t79 = t83 * t92 - t93;
	t78 = t83 * t93 - t92;
	t77 = t83 * t94 + t91;
	t1 = [t80, t89 * t82, t79, t80 * pkin(4) + t79 * qJ(5) + t85 * t87 + t90 * t89 + 0; t78, t87 * t82, t77, t78 * pkin(4) + t77 * qJ(5) - t89 * t85 + t90 * t87 + 0; t82 * t88, -t83, t82 * t86, sin(pkin(8)) * pkin(2) - t83 * pkin(7) + pkin(5) + 0 + (pkin(4) * t88 + qJ(5) * t86 + pkin(3)) * t82; 0, 0, 0, 1;];
	Tc_mdh = t1;
end