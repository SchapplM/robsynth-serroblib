% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:05
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:54
	% EndTime: 2020-11-04 22:05:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:54
	% EndTime: 2020-11-04 22:05:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t74 = cos(qJ(1));
	t73 = sin(qJ(1));
	t1 = [t74, -t73, 0, 0; t73, t74, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:54
	% EndTime: 2020-11-04 22:05:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t75 = sin(pkin(6));
	t78 = sin(qJ(1));
	t86 = t78 * t75;
	t77 = sin(qJ(2));
	t85 = t78 * t77;
	t79 = cos(qJ(2));
	t84 = t78 * t79;
	t80 = cos(qJ(1));
	t83 = t80 * t75;
	t82 = t80 * t77;
	t81 = t80 * t79;
	t76 = cos(pkin(6));
	t1 = [-t76 * t85 + t81, -t76 * t84 - t82, t86, t80 * pkin(1) + pkin(8) * t86 + 0; t76 * t82 + t84, t76 * t81 - t85, -t83, t78 * pkin(1) - pkin(8) * t83 + 0; t75 * t77, t75 * t79, t76, t76 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:54
	% EndTime: 2020-11-04 22:05:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->18), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t91 = sin(qJ(2));
	t92 = sin(qJ(1));
	t99 = t92 * t91;
	t93 = cos(qJ(2));
	t98 = t92 * t93;
	t94 = cos(qJ(1));
	t97 = t94 * t91;
	t96 = t94 * t93;
	t95 = pkin(2) * t91 - qJ(3) * t93;
	t90 = cos(pkin(6));
	t89 = sin(pkin(6));
	t88 = t93 * pkin(2) + t91 * qJ(3) + pkin(1);
	t87 = t89 * pkin(8) - t95 * t90;
	t1 = [t92 * t89, t90 * t99 - t96, t90 * t98 + t97, t87 * t92 + t88 * t94 + 0; -t94 * t89, -t90 * t97 - t98, -t90 * t96 + t99, -t87 * t94 + t88 * t92 + 0; t90, -t89 * t91, -t89 * t93, t90 * pkin(8) + t95 * t89 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:54
	% EndTime: 2020-11-04 22:05:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (36->23), mult. (60->36), div. (0->0), fcn. (81->8), ass. (0->23)
	t105 = sin(pkin(6));
	t110 = sin(qJ(1));
	t121 = t105 * t110;
	t111 = cos(qJ(2));
	t120 = t105 * t111;
	t112 = cos(qJ(1));
	t119 = t105 * t112;
	t109 = sin(qJ(2));
	t118 = t110 * t109;
	t117 = t110 * t111;
	t116 = t112 * t109;
	t115 = t112 * t111;
	t108 = pkin(2) + qJ(4);
	t114 = qJ(3) * t111 - t108 * t109;
	t113 = pkin(3) + pkin(8);
	t107 = cos(pkin(6));
	t106 = cos(pkin(11));
	t104 = sin(pkin(11));
	t103 = t109 * qJ(3) + t108 * t111 + pkin(1);
	t102 = t107 * t117 + t116;
	t101 = -t107 * t115 + t118;
	t100 = t105 * t113 + t114 * t107;
	t1 = [t102 * t104 + t106 * t121, t102 * t106 - t104 * t121, -t107 * t118 + t115, t100 * t110 + t103 * t112 + 0; t101 * t104 - t106 * t119, t101 * t106 + t104 * t119, t107 * t116 + t117, -t100 * t112 + t103 * t110 + 0; -t104 * t120 + t107 * t106, -t107 * t104 - t106 * t120, t105 * t109, -t114 * t105 + t113 * t107 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:54
	% EndTime: 2020-11-04 22:05:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (61->26), mult. (68->38), div. (0->0), fcn. (89->10), ass. (0->25)
	t132 = sin(pkin(6));
	t135 = sin(qJ(1));
	t145 = t132 * t135;
	t136 = cos(qJ(2));
	t144 = t132 * t136;
	t137 = cos(qJ(1));
	t143 = t132 * t137;
	t134 = sin(qJ(2));
	t142 = t135 * t134;
	t141 = t135 * t136;
	t140 = t137 * t134;
	t139 = t137 * t136;
	t127 = sin(pkin(11)) * pkin(4) + qJ(3);
	t130 = qJ(4) + pkin(2) + pkin(9);
	t138 = t127 * t136 - t130 * t134;
	t133 = cos(pkin(6));
	t131 = pkin(11) + qJ(5);
	t129 = cos(t131);
	t128 = sin(t131);
	t126 = cos(pkin(11)) * pkin(4) + pkin(3) + pkin(8);
	t125 = t133 * t141 + t140;
	t124 = -t133 * t139 + t142;
	t123 = t127 * t134 + t130 * t136 + pkin(1);
	t122 = t132 * t126 + t138 * t133;
	t1 = [t125 * t128 + t129 * t145, t125 * t129 - t128 * t145, -t133 * t142 + t139, t122 * t135 + t123 * t137 + 0; t124 * t128 - t129 * t143, t124 * t129 + t128 * t143, t133 * t140 + t141, -t122 * t137 + t123 * t135 + 0; -t128 * t144 + t133 * t129, -t133 * t128 - t129 * t144, t132 * t134, t126 * t133 - t138 * t132 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:54
	% EndTime: 2020-11-04 22:05:54
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (108->39), mult. (138->56), div. (0->0), fcn. (182->12), ass. (0->36)
	t164 = sin(pkin(6));
	t167 = sin(qJ(2));
	t180 = t164 * t167;
	t168 = sin(qJ(1));
	t179 = t164 * t168;
	t170 = cos(qJ(2));
	t178 = t164 * t170;
	t171 = cos(qJ(1));
	t177 = t164 * t171;
	t176 = t168 * t167;
	t175 = t168 * t170;
	t174 = t171 * t167;
	t173 = t171 * t170;
	t159 = sin(pkin(11)) * pkin(4) + qJ(3);
	t162 = qJ(4) + pkin(2) + pkin(9);
	t172 = t159 * t170 - t162 * t167;
	t169 = cos(qJ(6));
	t166 = sin(qJ(6));
	t165 = cos(pkin(6));
	t163 = pkin(11) + qJ(5);
	t161 = cos(t163);
	t160 = sin(t163);
	t158 = cos(pkin(11)) * pkin(4) + pkin(3) + pkin(8);
	t157 = -t165 * t176 + t173;
	t156 = t165 * t175 + t174;
	t155 = t165 * t174 + t175;
	t154 = -t165 * t173 + t176;
	t153 = t159 * t167 + t162 * t170 + pkin(1);
	t152 = -t160 * t178 + t161 * t165;
	t151 = t160 * t165 + t161 * t178;
	t150 = t154 * t160 - t161 * t177;
	t149 = t154 * t161 + t160 * t177;
	t148 = t156 * t160 + t161 * t179;
	t147 = -t156 * t161 + t160 * t179;
	t146 = t164 * t158 + t165 * t172;
	t1 = [t148 * t169 + t157 * t166, -t148 * t166 + t157 * t169, t147, pkin(5) * t148 + pkin(10) * t147 + t146 * t168 + t153 * t171 + 0; t150 * t169 + t155 * t166, -t150 * t166 + t155 * t169, -t149, pkin(5) * t150 - pkin(10) * t149 - t146 * t171 + t153 * t168 + 0; t152 * t169 + t166 * t180, -t152 * t166 + t169 * t180, t151, t152 * pkin(5) + t151 * pkin(10) + t158 * t165 - t164 * t172 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end