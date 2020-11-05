% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:11
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:11:17
	% EndTime: 2020-11-04 22:11:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:11:17
	% EndTime: 2020-11-04 22:11:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t82 = cos(qJ(1));
	t81 = sin(qJ(1));
	t1 = [t82, -t81, 0, 0; t81, t82, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:11:17
	% EndTime: 2020-11-04 22:11:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t83 = sin(pkin(6));
	t86 = sin(qJ(1));
	t94 = t86 * t83;
	t85 = sin(qJ(2));
	t93 = t86 * t85;
	t87 = cos(qJ(2));
	t92 = t86 * t87;
	t88 = cos(qJ(1));
	t91 = t88 * t83;
	t90 = t88 * t85;
	t89 = t88 * t87;
	t84 = cos(pkin(6));
	t1 = [-t84 * t93 + t89, -t84 * t92 - t90, t94, t88 * pkin(1) + pkin(8) * t94 + 0; t84 * t90 + t92, t84 * t89 - t93, -t91, t86 * pkin(1) - pkin(8) * t91 + 0; t83 * t85, t83 * t87, t84, t84 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:11:17
	% EndTime: 2020-11-04 22:11:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->22), mult. (60->36), div. (0->0), fcn. (81->8), ass. (0->21)
	t100 = sin(pkin(6));
	t103 = sin(qJ(2));
	t114 = t100 * t103;
	t104 = sin(qJ(1));
	t113 = t100 * t104;
	t106 = cos(qJ(1));
	t112 = t100 * t106;
	t111 = t104 * t103;
	t105 = cos(qJ(2));
	t110 = t104 * t105;
	t109 = t106 * t103;
	t108 = t106 * t105;
	t107 = pkin(2) * t103 - qJ(3) * t105;
	t102 = cos(pkin(6));
	t101 = cos(pkin(11));
	t99 = sin(pkin(11));
	t98 = t105 * pkin(2) + t103 * qJ(3) + pkin(1);
	t97 = -t102 * t111 + t108;
	t96 = t102 * t109 + t110;
	t95 = t100 * pkin(8) - t107 * t102;
	t1 = [t97 * t101 + t99 * t113, t101 * t113 - t97 * t99, t102 * t110 + t109, t95 * t104 + t98 * t106 + 0; t96 * t101 - t99 * t112, -t101 * t112 - t96 * t99, -t102 * t108 + t111, t98 * t104 - t95 * t106 + 0; t101 * t114 + t102 * t99, t102 * t101 - t99 * t114, -t100 * t105, t102 * pkin(8) + t107 * t100 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:11:17
	% EndTime: 2020-11-04 22:11:17
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (54->26), mult. (68->38), div. (0->0), fcn. (89->10), ass. (0->25)
	t124 = sin(pkin(6));
	t127 = sin(qJ(2));
	t138 = t124 * t127;
	t128 = sin(qJ(1));
	t137 = t124 * t128;
	t130 = cos(qJ(1));
	t136 = t124 * t130;
	t135 = t128 * t127;
	t129 = cos(qJ(2));
	t134 = t128 * t129;
	t133 = t130 * t127;
	t132 = t130 * t129;
	t120 = cos(pkin(11)) * pkin(3) + pkin(2);
	t126 = qJ(3) + pkin(9);
	t131 = t120 * t127 - t126 * t129;
	t125 = cos(pkin(6));
	t123 = pkin(11) + qJ(4);
	t122 = cos(t123);
	t121 = sin(t123);
	t119 = sin(pkin(11)) * pkin(3) + pkin(8);
	t118 = -t125 * t135 + t132;
	t117 = t125 * t133 + t134;
	t116 = t120 * t129 + t126 * t127 + pkin(1);
	t115 = t124 * t119 - t131 * t125;
	t1 = [t118 * t122 + t121 * t137, -t118 * t121 + t122 * t137, t125 * t134 + t133, t115 * t128 + t116 * t130 + 0; t117 * t122 - t121 * t136, -t117 * t121 - t122 * t136, -t125 * t132 + t135, -t115 * t130 + t116 * t128 + 0; t125 * t121 + t122 * t138, -t121 * t138 + t125 * t122, -t124 * t129, t119 * t125 + t131 * t124 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:11:17
	% EndTime: 2020-11-04 22:11:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (82->36), mult. (104->44), div. (0->0), fcn. (135->10), ass. (0->31)
	t154 = sin(pkin(6));
	t157 = sin(qJ(2));
	t168 = t154 * t157;
	t158 = sin(qJ(1));
	t167 = t154 * t158;
	t160 = cos(qJ(1));
	t166 = t154 * t160;
	t165 = t158 * t157;
	t159 = cos(qJ(2));
	t164 = t158 * t159;
	t163 = t160 * t157;
	t162 = t160 * t159;
	t150 = cos(pkin(11)) * pkin(3) + pkin(2);
	t156 = qJ(3) + pkin(9);
	t161 = t150 * t157 - t156 * t159;
	t155 = cos(pkin(6));
	t153 = pkin(11) + qJ(4);
	t152 = cos(t153);
	t151 = sin(t153);
	t149 = sin(pkin(11)) * pkin(3) + pkin(8);
	t148 = -t155 * t165 + t162;
	t147 = t155 * t163 + t164;
	t146 = t150 * t159 + t156 * t157 + pkin(1);
	t145 = t155 * t151 + t152 * t168;
	t144 = t151 * t168 - t155 * t152;
	t143 = t154 * t149 - t161 * t155;
	t142 = -t148 * t151 + t152 * t167;
	t141 = t148 * t152 + t151 * t167;
	t140 = t147 * t152 - t151 * t166;
	t139 = t147 * t151 + t152 * t166;
	t1 = [t155 * t164 + t163, -t141, -t142, t141 * pkin(4) - t142 * qJ(5) + t143 * t158 + t146 * t160 + 0; -t155 * t162 + t165, -t140, t139, t140 * pkin(4) + t139 * qJ(5) - t143 * t160 + t146 * t158 + 0; -t154 * t159, -t145, t144, t145 * pkin(4) + t144 * qJ(5) + t149 * t155 + t161 * t154 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:11:17
	% EndTime: 2020-11-04 22:11:17
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (125->38), mult. (156->58), div. (0->0), fcn. (190->14), ass. (0->40)
	t183 = sin(pkin(11));
	t185 = cos(pkin(11));
	t195 = pkin(4) + pkin(10);
	t177 = qJ(5) * t183 + t195 * t185;
	t178 = qJ(5) * t185 - t183 * t195;
	t188 = sin(qJ(4));
	t192 = cos(qJ(4));
	t209 = -t183 * pkin(3) - t177 * t188 + t178 * t192 - pkin(8);
	t184 = sin(pkin(6));
	t189 = sin(qJ(2));
	t207 = t184 * t189;
	t190 = sin(qJ(1));
	t206 = t184 * t190;
	t193 = cos(qJ(2));
	t205 = t184 * t193;
	t194 = cos(qJ(1));
	t204 = t184 * t194;
	t203 = t190 * t189;
	t202 = t190 * t193;
	t201 = t194 * t189;
	t200 = t194 * t193;
	t171 = t185 * pkin(3) + t177 * t192 + t178 * t188 + pkin(2);
	t181 = qJ(3) + pkin(5) + pkin(9);
	t199 = t171 * t189 - t181 * t193;
	t186 = cos(pkin(6));
	t174 = t186 * t201 + t202;
	t182 = pkin(11) + qJ(4);
	t179 = sin(t182);
	t180 = cos(t182);
	t197 = t174 * t179 + t180 * t204;
	t176 = -t186 * t203 + t200;
	t196 = t176 * t179 - t180 * t206;
	t191 = cos(qJ(6));
	t187 = sin(qJ(6));
	t175 = t186 * t202 + t201;
	t173 = t186 * t200 - t203;
	t172 = t179 * t207 - t186 * t180;
	t170 = t171 * t193 + t181 * t189 + pkin(1);
	t169 = -t184 * t209 - t199 * t186;
	t1 = [t175 * t191 + t196 * t187, -t175 * t187 + t196 * t191, t176 * t180 + t179 * t206, t169 * t190 + t170 * t194 + 0; -t191 * t173 + t197 * t187, t187 * t173 + t197 * t191, t174 * t180 - t179 * t204, -t169 * t194 + t170 * t190 + 0; t172 * t187 - t191 * t205, t172 * t191 + t187 * t205, t186 * t179 + t180 * t207, t199 * t184 - t209 * t186 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end