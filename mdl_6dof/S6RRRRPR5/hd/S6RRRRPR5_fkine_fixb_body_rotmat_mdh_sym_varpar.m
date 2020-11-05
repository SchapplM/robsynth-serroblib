% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:38
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:38:41
	% EndTime: 2020-11-04 22:38:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:38:41
	% EndTime: 2020-11-04 22:38:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t70 = cos(qJ(1));
	t69 = sin(qJ(1));
	t1 = [t70, -t69, 0, 0; t69, t70, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:38:41
	% EndTime: 2020-11-04 22:38:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t74 = cos(qJ(1));
	t73 = cos(qJ(2));
	t72 = sin(qJ(1));
	t71 = sin(qJ(2));
	t1 = [t74 * t73, -t74 * t71, t72, t74 * pkin(1) + t72 * pkin(7) + 0; t72 * t73, -t72 * t71, -t74, t72 * pkin(1) - t74 * pkin(7) + 0; t71, t73, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:38:41
	% EndTime: 2020-11-04 22:38:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t81 = pkin(8) + pkin(7);
	t80 = cos(qJ(1));
	t79 = sin(qJ(1));
	t78 = qJ(2) + qJ(3);
	t77 = cos(t78);
	t76 = sin(t78);
	t75 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t80 * t77, -t80 * t76, t79, t80 * t75 + t81 * t79 + 0; t79 * t77, -t79 * t76, -t80, t79 * t75 - t80 * t81 + 0; t76, t77, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:38:41
	% EndTime: 2020-11-04 22:38:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t86 = sin(qJ(4));
	t87 = sin(qJ(1));
	t95 = t87 * t86;
	t88 = cos(qJ(4));
	t94 = t87 * t88;
	t89 = cos(qJ(1));
	t93 = t89 * t86;
	t92 = t89 * t88;
	t85 = qJ(2) + qJ(3);
	t83 = sin(t85);
	t84 = cos(t85);
	t91 = pkin(3) * t84 + pkin(9) * t83 + cos(qJ(2)) * pkin(2) + pkin(1);
	t90 = pkin(8) + pkin(7);
	t1 = [t84 * t92 + t95, -t84 * t93 + t94, t89 * t83, t90 * t87 + t91 * t89 + 0; t84 * t94 - t93, -t84 * t95 - t92, t87 * t83, t91 * t87 - t89 * t90 + 0; t83 * t88, -t83 * t86, -t84, t83 * pkin(3) - t84 * pkin(9) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:38:41
	% EndTime: 2020-11-04 22:38:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (52->24), mult. (53->28), div. (0->0), fcn. (70->8), ass. (0->18)
	t104 = sin(qJ(4));
	t105 = sin(qJ(1));
	t113 = t105 * t104;
	t106 = cos(qJ(4));
	t112 = t105 * t106;
	t107 = cos(qJ(1));
	t111 = t107 * t104;
	t110 = t107 * t106;
	t103 = qJ(2) + qJ(3);
	t101 = sin(t103);
	t102 = cos(t103);
	t109 = pkin(3) * t102 + pkin(9) * t101 + cos(qJ(2)) * pkin(2) + pkin(1);
	t108 = pkin(8) + pkin(7);
	t99 = t102 * t110 + t113;
	t98 = t102 * t111 - t112;
	t97 = t102 * t112 - t111;
	t96 = t102 * t113 + t110;
	t1 = [t99, t107 * t101, t98, t99 * pkin(4) + t98 * qJ(5) + t108 * t105 + t109 * t107 + 0; t97, t105 * t101, t96, t97 * pkin(4) + t96 * qJ(5) + t109 * t105 - t107 * t108 + 0; t101 * t106, -t102, t101 * t104, sin(qJ(2)) * pkin(2) - t102 * pkin(9) + pkin(6) + 0 + (pkin(4) * t106 + qJ(5) * t104 + pkin(3)) * t101; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:38:41
	% EndTime: 2020-11-04 22:38:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (84->35), mult. (81->39), div. (0->0), fcn. (97->16), ass. (0->26)
	t123 = sin(qJ(6));
	t124 = sin(qJ(4));
	t128 = cos(qJ(6));
	t129 = cos(qJ(4));
	t115 = -t123 * t129 + t128 * t124;
	t127 = sin(qJ(1));
	t138 = t127 * t115;
	t116 = t123 * t124 + t128 * t129;
	t137 = t127 * t116;
	t131 = cos(qJ(1));
	t136 = t131 * t115;
	t135 = t131 * t116;
	t122 = qJ(2) + qJ(3);
	t133 = pkin(4) + pkin(5);
	t134 = qJ(5) * t129 - t133 * t124 - pkin(7) - pkin(8);
	t132 = pkin(10) - pkin(9);
	t130 = cos(qJ(3));
	t126 = sin(qJ(2));
	t125 = sin(qJ(3));
	t121 = -qJ(4) + t122;
	t120 = qJ(4) + t122;
	t119 = cos(t122);
	t118 = sin(t122);
	t117 = qJ(5) * t124 + t133 * t129 + pkin(3);
	t114 = (t117 * t130 - t132 * t125 + pkin(2)) * cos(qJ(2)) + pkin(1) + (-t117 * t125 - t132 * t130) * t126;
	t1 = [t119 * t135 + t138, t119 * t136 - t137, -t131 * t118, t114 * t131 - t134 * t127 + 0; t119 * t137 - t136, t119 * t138 + t135, -t127 * t118, t114 * t127 + t134 * t131 + 0; t118 * t116, t118 * t115, t119, t132 * t119 + t126 * pkin(2) + t118 * pkin(3) + 0 + pkin(6) + (sin(t121) / 0.2e1 + sin(t120) / 0.2e1) * t133 + (cos(t121) / 0.2e1 - cos(t120) / 0.2e1) * qJ(5); 0, 0, 0, 1;];
	Tc_mdh = t1;
end