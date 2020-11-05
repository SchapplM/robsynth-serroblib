% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR16 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:39
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRR16_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR16_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:22
	% EndTime: 2020-11-04 20:39:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:22
	% EndTime: 2020-11-04 20:39:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t71 = cos(qJ(1));
	t70 = sin(qJ(1));
	t1 = [t71, -t70, 0, 0; t70, t71, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:22
	% EndTime: 2020-11-04 20:39:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t72 = sin(pkin(5));
	t75 = sin(qJ(1));
	t83 = t75 * t72;
	t74 = sin(qJ(2));
	t82 = t75 * t74;
	t76 = cos(qJ(2));
	t81 = t75 * t76;
	t77 = cos(qJ(1));
	t80 = t77 * t72;
	t79 = t77 * t74;
	t78 = t77 * t76;
	t73 = cos(pkin(5));
	t1 = [-t73 * t82 + t78, -t73 * t81 - t79, t83, t77 * pkin(1) + pkin(7) * t83 + 0; t73 * t79 + t81, t73 * t78 - t82, -t80, t75 * pkin(1) - pkin(7) * t80 + 0; t72 * t74, t72 * t76, t73, t73 * pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:22
	% EndTime: 2020-11-04 20:39:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->18), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t88 = sin(qJ(2));
	t89 = sin(qJ(1));
	t96 = t89 * t88;
	t90 = cos(qJ(2));
	t95 = t89 * t90;
	t91 = cos(qJ(1));
	t94 = t91 * t88;
	t93 = t91 * t90;
	t92 = pkin(2) * t88 - qJ(3) * t90;
	t87 = cos(pkin(5));
	t86 = sin(pkin(5));
	t85 = t90 * pkin(2) + t88 * qJ(3) + pkin(1);
	t84 = t86 * pkin(7) - t92 * t87;
	t1 = [t89 * t86, t87 * t96 - t93, t87 * t95 + t94, t84 * t89 + t85 * t91 + 0; -t91 * t86, -t87 * t94 - t95, -t87 * t93 + t96, -t84 * t91 + t85 * t89 + 0; t87, -t86 * t88, -t86 * t90, t87 * pkin(7) + t92 * t86 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:22
	% EndTime: 2020-11-04 20:39:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->24), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->22)
	t100 = sin(pkin(5));
	t105 = cos(qJ(4));
	t117 = t100 * t105;
	t102 = sin(qJ(4));
	t116 = t102 * t100;
	t103 = sin(qJ(2));
	t104 = sin(qJ(1));
	t115 = t104 * t103;
	t106 = cos(qJ(2));
	t114 = t104 * t106;
	t113 = t105 * t106;
	t107 = cos(qJ(1));
	t112 = t107 * t103;
	t111 = t107 * t106;
	t109 = pkin(2) + pkin(8);
	t110 = qJ(3) * t106 - t103 * t109;
	t108 = pkin(3) + pkin(7);
	t101 = cos(pkin(5));
	t99 = t103 * qJ(3) + t109 * t106 + pkin(1);
	t98 = -t101 * t111 + t115;
	t97 = t100 * t108 + t110 * t101;
	t1 = [t104 * t117 + (t101 * t114 + t112) * t102, (t101 * t113 - t116) * t104 + t105 * t112, -t101 * t115 + t111, t97 * t104 + t99 * t107 + 0; t98 * t102 - t107 * t117, t98 * t105 + t107 * t116, t101 * t112 + t114, t99 * t104 - t97 * t107 + 0; t101 * t105 - t106 * t116, -t100 * t113 - t101 * t102, t100 * t103, -t110 * t100 + t108 * t101 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:22
	% EndTime: 2020-11-04 20:39:22
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (65->36), mult. (112->59), div. (0->0), fcn. (146->10), ass. (0->29)
	t126 = sin(qJ(5));
	t128 = sin(qJ(2));
	t145 = t126 * t128;
	t124 = sin(pkin(5));
	t127 = sin(qJ(4));
	t144 = t127 * t124;
	t132 = cos(qJ(2));
	t143 = t127 * t132;
	t130 = cos(qJ(5));
	t142 = t128 * t130;
	t133 = cos(qJ(1));
	t141 = t128 * t133;
	t129 = sin(qJ(1));
	t140 = t129 * t128;
	t131 = cos(qJ(4));
	t139 = t131 * t132;
	t138 = t133 * t132;
	t123 = t127 * pkin(4) - pkin(9) * t131 + qJ(3);
	t134 = pkin(2) + pkin(8);
	t137 = t123 * t132 - t134 * t128;
	t136 = t131 * pkin(4) + t127 * pkin(9) + pkin(3) + pkin(7);
	t125 = cos(pkin(5));
	t120 = t131 * t124 + t125 * t143;
	t135 = t120 * t126 + t125 * t142;
	t122 = -t127 * t145 + t130 * t132;
	t121 = -t124 * t143 + t125 * t131;
	t119 = t123 * t128 + t134 * t132 + pkin(1);
	t118 = t124 * t136 + t137 * t125;
	t1 = [(t120 * t129 + t127 * t141) * t130 + (-t125 * t140 + t138) * t126, t133 * t122 - t135 * t129, (-t125 * t139 + t144) * t129 - t131 * t141, t118 * t129 + t119 * t133 + 0; (-t120 * t130 + t125 * t145) * t133 + t129 * (t126 * t132 + t127 * t142), t129 * t122 + t135 * t133, (t125 * t138 - t140) * t131 - t133 * t144, -t118 * t133 + t119 * t129 + 0; t121 * t130 + t124 * t145, -t121 * t126 + t124 * t142, t124 * t139 + t125 * t127, -t137 * t124 + t136 * t125 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end