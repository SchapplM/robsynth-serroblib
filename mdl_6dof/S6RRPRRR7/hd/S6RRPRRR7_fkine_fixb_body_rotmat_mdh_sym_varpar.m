% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:18
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:58
	% EndTime: 2020-11-04 22:18:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:58
	% EndTime: 2020-11-04 22:18:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t66 = cos(qJ(1));
	t65 = sin(qJ(1));
	t1 = [t66, -t65, 0, 0; t65, t66, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:58
	% EndTime: 2020-11-04 22:18:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t70 = cos(qJ(1));
	t69 = cos(qJ(2));
	t68 = sin(qJ(1));
	t67 = sin(qJ(2));
	t1 = [t70 * t69, -t70 * t67, t68, t70 * pkin(1) + t68 * pkin(7) + 0; t68 * t69, -t68 * t67, -t70, t68 * pkin(1) - t70 * pkin(7) + 0; t67, t69, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:58
	% EndTime: 2020-11-04 22:18:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t75 = cos(qJ(1));
	t74 = cos(qJ(2));
	t73 = sin(qJ(1));
	t72 = sin(qJ(2));
	t71 = t74 * pkin(2) + t72 * qJ(3) + pkin(1);
	t1 = [t75 * t74, t73, t75 * t72, t73 * pkin(7) + t71 * t75 + 0; t73 * t74, -t75, t73 * t72, -t75 * pkin(7) + t71 * t73 + 0; t72, 0, -t74, t72 * pkin(2) - t74 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:58
	% EndTime: 2020-11-04 22:18:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->16), mult. (26->16), div. (0->0), fcn. (40->6), ass. (0->12)
	t88 = cos(qJ(4));
	t87 = sin(qJ(4));
	t86 = pkin(2) + pkin(3);
	t85 = pkin(7) - pkin(8);
	t84 = cos(qJ(1));
	t83 = cos(qJ(2));
	t82 = sin(qJ(1));
	t81 = sin(qJ(2));
	t78 = t81 * qJ(3) + t86 * t83 + pkin(1);
	t77 = t81 * t88 - t83 * t87;
	t76 = -t81 * t87 - t83 * t88;
	t1 = [-t84 * t76, t84 * t77, -t82, t78 * t84 + t85 * t82 + 0; -t82 * t76, t82 * t77, t84, t78 * t82 - t85 * t84 + 0; t77, t76, 0, -t83 * qJ(3) + t86 * t81 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:58
	% EndTime: 2020-11-04 22:18:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (41->24), mult. (56->28), div. (0->0), fcn. (78->8), ass. (0->19)
	t95 = sin(qJ(5));
	t98 = sin(qJ(1));
	t107 = t98 * t95;
	t99 = cos(qJ(5));
	t106 = t98 * t99;
	t102 = cos(qJ(1));
	t105 = t102 * t95;
	t104 = t102 * t99;
	t103 = pkin(7) - pkin(8);
	t101 = cos(qJ(2));
	t100 = cos(qJ(4));
	t97 = sin(qJ(2));
	t96 = sin(qJ(4));
	t93 = -t96 * pkin(4) + t100 * pkin(9) - qJ(3);
	t92 = t100 * pkin(4) + t96 * pkin(9) + pkin(2) + pkin(3);
	t91 = t101 * t100 + t97 * t96;
	t90 = t97 * t100 - t101 * t96;
	t89 = t92 * t101 - t93 * t97 + pkin(1);
	t1 = [t91 * t104 - t107, -t91 * t105 - t106, -t102 * t90, t89 * t102 + t103 * t98 + 0; t91 * t106 + t105, -t91 * t107 + t104, -t98 * t90, -t103 * t102 + t89 * t98 + 0; t90 * t99, -t90 * t95, t91, t93 * t101 + t92 * t97 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:58
	% EndTime: 2020-11-04 22:18:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (65->28), mult. (64->30), div. (0->0), fcn. (86->10), ass. (0->22)
	t117 = qJ(5) + qJ(6);
	t115 = sin(t117);
	t120 = sin(qJ(1));
	t129 = t120 * t115;
	t116 = cos(t117);
	t128 = t120 * t116;
	t123 = cos(qJ(1));
	t127 = t123 * t115;
	t126 = t123 * t116;
	t114 = cos(qJ(5)) * pkin(5) + pkin(4);
	t118 = sin(qJ(4));
	t121 = cos(qJ(4));
	t124 = pkin(10) + pkin(9);
	t125 = t114 * t118 - t124 * t121 + qJ(3);
	t122 = cos(qJ(2));
	t119 = sin(qJ(2));
	t112 = sin(qJ(5)) * pkin(5) - pkin(7) + pkin(8);
	t111 = t119 * t118 + t122 * t121;
	t110 = -t122 * t118 + t119 * t121;
	t109 = t114 * t121 + t124 * t118 + pkin(2) + pkin(3);
	t108 = t109 * t122 + t125 * t119 + pkin(1);
	t1 = [t111 * t126 - t129, -t111 * t127 - t128, -t123 * t110, t108 * t123 - t112 * t120 + 0; t111 * t128 + t127, -t111 * t129 + t126, -t120 * t110, t108 * t120 + t112 * t123 + 0; t110 * t116, -t110 * t115, t111, t109 * t119 - t125 * t122 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end