% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:10
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:32
	% EndTime: 2020-11-04 22:10:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:32
	% EndTime: 2020-11-04 22:10:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t69 = cos(qJ(1));
	t68 = sin(qJ(1));
	t1 = [t69, -t68, 0, 0; t68, t69, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:32
	% EndTime: 2020-11-04 22:10:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t73 = cos(qJ(1));
	t72 = cos(qJ(2));
	t71 = sin(qJ(1));
	t70 = sin(qJ(2));
	t1 = [t73 * t72, -t73 * t70, t71, t73 * pkin(1) + t71 * pkin(7) + 0; t71 * t72, -t71 * t70, -t73, t71 * pkin(1) - t73 * pkin(7) + 0; t70, t72, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:32
	% EndTime: 2020-11-04 22:10:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t78 = sin(qJ(1));
	t79 = cos(qJ(2));
	t83 = t78 * t79;
	t75 = sin(pkin(10));
	t80 = cos(qJ(1));
	t82 = t80 * t75;
	t76 = cos(pkin(10));
	t81 = t80 * t76;
	t77 = sin(qJ(2));
	t74 = t79 * pkin(2) + t77 * qJ(3) + pkin(1);
	t1 = [t78 * t75 + t79 * t81, t78 * t76 - t79 * t82, t80 * t77, t78 * pkin(7) + t74 * t80 + 0; t76 * t83 - t82, -t75 * t83 - t81, t78 * t77, -t80 * pkin(7) + t74 * t78 + 0; t77 * t76, -t77 * t75, -t79, t77 * pkin(2) - t79 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:32
	% EndTime: 2020-11-04 22:10:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t92 = sin(qJ(1));
	t93 = cos(qJ(2));
	t97 = t92 * t93;
	t89 = pkin(10) + qJ(4);
	t87 = sin(t89);
	t94 = cos(qJ(1));
	t96 = t94 * t87;
	t88 = cos(t89);
	t95 = t94 * t88;
	t91 = sin(qJ(2));
	t90 = qJ(3) + pkin(8);
	t86 = cos(pkin(10)) * pkin(3) + pkin(2);
	t85 = sin(pkin(10)) * pkin(3) + pkin(7);
	t84 = t86 * t93 + t90 * t91 + pkin(1);
	t1 = [t92 * t87 + t93 * t95, t92 * t88 - t93 * t96, t94 * t91, t84 * t94 + t85 * t92 + 0; t88 * t97 - t96, -t87 * t97 - t95, t92 * t91, t84 * t92 - t85 * t94 + 0; t91 * t88, -t91 * t87, -t93, t91 * t86 - t93 * t90 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:32
	% EndTime: 2020-11-04 22:10:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (54->24), mult. (51->29), div. (0->0), fcn. (68->8), ass. (0->19)
	t110 = sin(qJ(1));
	t111 = cos(qJ(2));
	t115 = t110 * t111;
	t107 = pkin(10) + qJ(4);
	t105 = sin(t107);
	t112 = cos(qJ(1));
	t114 = t112 * t105;
	t106 = cos(t107);
	t113 = t112 * t106;
	t109 = sin(qJ(2));
	t108 = qJ(3) + pkin(8);
	t104 = cos(pkin(10)) * pkin(3) + pkin(2);
	t103 = sin(pkin(10)) * pkin(3) + pkin(7);
	t102 = t104 * t111 + t108 * t109 + pkin(1);
	t101 = t110 * t105 + t111 * t113;
	t100 = -t110 * t106 + t111 * t114;
	t99 = t106 * t115 - t114;
	t98 = t105 * t115 + t113;
	t1 = [t101, t112 * t109, t100, t101 * pkin(4) + t100 * qJ(5) + t102 * t112 + t103 * t110 + 0; t99, t110 * t109, t98, t99 * pkin(4) + t98 * qJ(5) + t102 * t110 - t103 * t112 + 0; t109 * t106, -t111, t109 * t105, -t111 * t108 + pkin(6) + 0 + (pkin(4) * t106 + qJ(5) * t105 + t104) * t109; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:32
	% EndTime: 2020-11-04 22:10:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (80->32), mult. (85->42), div. (0->0), fcn. (108->12), ass. (0->29)
	t128 = sin(pkin(10));
	t129 = cos(pkin(10));
	t138 = pkin(4) + pkin(5);
	t122 = qJ(5) * t128 + t138 * t129;
	t123 = qJ(5) * t129 - t128 * t138;
	t131 = sin(qJ(4));
	t135 = cos(qJ(4));
	t145 = -t128 * pkin(3) - t122 * t131 + t123 * t135 - pkin(7);
	t130 = sin(qJ(6));
	t137 = cos(qJ(1));
	t143 = t130 * t137;
	t133 = sin(qJ(1));
	t142 = t133 * t130;
	t134 = cos(qJ(6));
	t141 = t134 * t133;
	t140 = t137 * t134;
	t136 = cos(qJ(2));
	t132 = sin(qJ(2));
	t127 = pkin(10) + qJ(4);
	t126 = qJ(3) + pkin(8) - pkin(9);
	t125 = cos(t127);
	t124 = sin(t127);
	t121 = t136 * t140 - t142;
	t120 = t136 * t141 + t143;
	t119 = t136 * t142 - t140;
	t118 = t136 * t143 + t141;
	t117 = t129 * pkin(3) + t122 * t135 + t123 * t131 + pkin(2);
	t116 = t117 * t136 + t126 * t132 + pkin(1);
	t1 = [t124 * t118 + t125 * t121, -t125 * t118 + t121 * t124, -t137 * t132, t116 * t137 - t145 * t133 + 0; t124 * t119 + t120 * t125, -t125 * t119 + t124 * t120, -t133 * t132, t116 * t133 + t145 * t137 + 0; t132 * (t130 * t124 + t125 * t134), t132 * (t124 * t134 - t130 * t125), t136, t117 * t132 - t126 * t136 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end