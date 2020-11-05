% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:01
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:29
	% EndTime: 2020-11-04 22:01:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:29
	% EndTime: 2020-11-04 22:01:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t72 = cos(qJ(1));
	t71 = sin(qJ(1));
	t1 = [t72, -t71, 0, 0; t71, t72, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:29
	% EndTime: 2020-11-04 22:01:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t76 = cos(qJ(1));
	t75 = cos(qJ(2));
	t74 = sin(qJ(1));
	t73 = sin(qJ(2));
	t1 = [t76 * t75, -t76 * t73, t74, t76 * pkin(1) + t74 * pkin(7) + 0; t74 * t75, -t74 * t73, -t76, t74 * pkin(1) - t76 * pkin(7) + 0; t73, t75, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:29
	% EndTime: 2020-11-04 22:01:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t81 = sin(qJ(1));
	t82 = cos(qJ(2));
	t86 = t81 * t82;
	t78 = sin(pkin(9));
	t83 = cos(qJ(1));
	t85 = t83 * t78;
	t79 = cos(pkin(9));
	t84 = t83 * t79;
	t80 = sin(qJ(2));
	t77 = t82 * pkin(2) + t80 * qJ(3) + pkin(1);
	t1 = [t81 * t78 + t82 * t84, t81 * t79 - t82 * t85, t83 * t80, t81 * pkin(7) + t77 * t83 + 0; t79 * t86 - t85, -t78 * t86 - t84, t81 * t80, -t83 * pkin(7) + t77 * t81 + 0; t80 * t79, -t80 * t78, -t82, t80 * pkin(2) - t82 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:30
	% EndTime: 2020-11-04 22:01:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->18), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->13)
	t93 = sin(qJ(1));
	t94 = cos(qJ(2));
	t98 = t93 * t94;
	t90 = sin(pkin(9));
	t95 = cos(qJ(1));
	t97 = t95 * t90;
	t91 = cos(pkin(9));
	t96 = t95 * t91;
	t92 = sin(qJ(2));
	t89 = -t90 * pkin(3) + qJ(4) * t91 - pkin(7);
	t88 = t91 * pkin(3) + t90 * qJ(4) + pkin(2);
	t87 = t92 * qJ(3) + t88 * t94 + pkin(1);
	t1 = [t93 * t90 + t94 * t96, t95 * t92, -t93 * t91 + t94 * t97, t87 * t95 - t89 * t93 + 0; t91 * t98 - t97, t93 * t92, t90 * t98 + t96, t87 * t93 + t89 * t95 + 0; t92 * t91, -t94, t92 * t90, -t94 * qJ(3) + t88 * t92 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:30
	% EndTime: 2020-11-04 22:01:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (45->23), mult. (56->30), div. (0->0), fcn. (79->8), ass. (0->18)
	t108 = sin(qJ(1));
	t110 = cos(qJ(2));
	t115 = t108 * t110;
	t111 = cos(qJ(1));
	t114 = t110 * t111;
	t103 = sin(pkin(9));
	t104 = cos(pkin(9));
	t112 = pkin(3) + pkin(4);
	t113 = qJ(4) * t104 - t112 * t103 - pkin(7);
	t109 = cos(qJ(5));
	t107 = sin(qJ(2));
	t106 = sin(qJ(5));
	t105 = qJ(3) - pkin(8);
	t102 = t103 * qJ(4) + t112 * t104 + pkin(2);
	t101 = t106 * t103 + t109 * t104;
	t100 = t109 * t103 - t106 * t104;
	t99 = t102 * t110 + t105 * t107 + pkin(1);
	t1 = [t108 * t100 + t101 * t114, t100 * t114 - t108 * t101, -t111 * t107, -t113 * t108 + t99 * t111 + 0; -t100 * t111 + t101 * t115, t100 * t115 + t101 * t111, -t108 * t107, t99 * t108 + t113 * t111 + 0; t107 * t101, t107 * t100, t110, t102 * t107 - t105 * t110 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:30
	% EndTime: 2020-11-04 22:01:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (65->30), mult. (86->38), div. (0->0), fcn. (109->8), ass. (0->20)
	t129 = sin(qJ(1));
	t131 = cos(qJ(2));
	t136 = t129 * t131;
	t132 = cos(qJ(1));
	t135 = t131 * t132;
	t124 = sin(pkin(9));
	t125 = cos(pkin(9));
	t120 = -pkin(5) * t125 + qJ(6) * t124;
	t121 = t124 * pkin(5) + qJ(6) * t125;
	t127 = sin(qJ(5));
	t130 = cos(qJ(5));
	t133 = pkin(3) + pkin(4);
	t134 = qJ(4) * t125 - t120 * t127 - t121 * t130 - t133 * t124 - pkin(7);
	t128 = sin(qJ(2));
	t126 = qJ(3) - pkin(8);
	t119 = t127 * t124 + t130 * t125;
	t118 = t130 * t124 - t127 * t125;
	t117 = t124 * qJ(4) - t120 * t130 + t121 * t127 + t133 * t125 + pkin(2);
	t116 = t117 * t131 + t126 * t128 + pkin(1);
	t1 = [t129 * t118 + t119 * t135, -t132 * t128, -t118 * t135 + t129 * t119, t116 * t132 - t134 * t129 + 0; -t118 * t132 + t119 * t136, -t129 * t128, -t118 * t136 - t119 * t132, t116 * t129 + t134 * t132 + 0; t128 * t119, t131, -t128 * t118, t117 * t128 - t126 * t131 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end