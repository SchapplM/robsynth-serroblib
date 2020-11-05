% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:20
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:03
	% EndTime: 2020-11-04 22:20:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:03
	% EndTime: 2020-11-04 22:20:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t78 = cos(qJ(1));
	t77 = sin(qJ(1));
	t1 = [t78, -t77, 0, 0; t77, t78, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:03
	% EndTime: 2020-11-04 22:20:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t79 = sin(pkin(6));
	t82 = sin(qJ(1));
	t90 = t82 * t79;
	t81 = sin(qJ(2));
	t89 = t82 * t81;
	t83 = cos(qJ(2));
	t88 = t82 * t83;
	t84 = cos(qJ(1));
	t87 = t84 * t79;
	t86 = t84 * t81;
	t85 = t84 * t83;
	t80 = cos(pkin(6));
	t1 = [-t80 * t89 + t85, -t80 * t88 - t86, t90, t84 * pkin(1) + pkin(8) * t90 + 0; t80 * t86 + t88, t80 * t85 - t89, -t87, t82 * pkin(1) - pkin(8) * t87 + 0; t79 * t81, t79 * t83, t80, t80 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:03
	% EndTime: 2020-11-04 22:20:03
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->22), mult. (60->36), div. (0->0), fcn. (81->8), ass. (0->21)
	t96 = sin(pkin(6));
	t99 = sin(qJ(2));
	t110 = t96 * t99;
	t100 = sin(qJ(1));
	t109 = t100 * t96;
	t108 = t100 * t99;
	t102 = cos(qJ(1));
	t107 = t102 * t96;
	t106 = t102 * t99;
	t101 = cos(qJ(2));
	t105 = t100 * t101;
	t104 = t102 * t101;
	t103 = pkin(2) * t99 - qJ(3) * t101;
	t98 = cos(pkin(6));
	t97 = cos(pkin(12));
	t95 = sin(pkin(12));
	t94 = t101 * pkin(2) + t99 * qJ(3) + pkin(1);
	t93 = t98 * t106 + t105;
	t92 = t98 * t108 - t104;
	t91 = t96 * pkin(8) - t103 * t98;
	t1 = [t95 * t109 - t92 * t97, t97 * t109 + t92 * t95, t98 * t105 + t106, t91 * t100 + t94 * t102 + 0; -t95 * t107 + t93 * t97, -t97 * t107 - t93 * t95, -t98 * t104 + t108, t94 * t100 - t91 * t102 + 0; t97 * t110 + t98 * t95, -t95 * t110 + t98 * t97, -t96 * t101, t98 * pkin(8) + t103 * t96 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:03
	% EndTime: 2020-11-04 22:20:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (54->26), mult. (68->38), div. (0->0), fcn. (89->10), ass. (0->25)
	t120 = sin(pkin(6));
	t123 = sin(qJ(2));
	t134 = t120 * t123;
	t124 = sin(qJ(1));
	t133 = t120 * t124;
	t126 = cos(qJ(1));
	t132 = t120 * t126;
	t131 = t124 * t123;
	t125 = cos(qJ(2));
	t130 = t124 * t125;
	t129 = t126 * t123;
	t128 = t126 * t125;
	t116 = cos(pkin(12)) * pkin(3) + pkin(2);
	t122 = qJ(3) + pkin(9);
	t127 = t116 * t123 - t122 * t125;
	t121 = cos(pkin(6));
	t119 = pkin(12) + qJ(4);
	t118 = cos(t119);
	t117 = sin(t119);
	t115 = sin(pkin(12)) * pkin(3) + pkin(8);
	t114 = t121 * t129 + t130;
	t113 = t121 * t131 - t128;
	t112 = t116 * t125 + t122 * t123 + pkin(1);
	t111 = t120 * t115 - t127 * t121;
	t1 = [-t113 * t118 + t117 * t133, t113 * t117 + t118 * t133, t121 * t130 + t129, t111 * t124 + t112 * t126 + 0; t114 * t118 - t117 * t132, -t114 * t117 - t118 * t132, -t121 * t128 + t131, -t111 * t126 + t112 * t124 + 0; t121 * t117 + t118 * t134, -t117 * t134 + t121 * t118, -t120 * t125, t115 * t121 + t127 * t120 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:03
	% EndTime: 2020-11-04 22:20:03
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (100->38), mult. (138->56), div. (0->0), fcn. (182->12), ass. (0->36)
	t152 = sin(pkin(6));
	t156 = sin(qJ(2));
	t169 = t152 * t156;
	t157 = sin(qJ(1));
	t168 = t152 * t157;
	t159 = cos(qJ(2));
	t167 = t152 * t159;
	t160 = cos(qJ(1));
	t166 = t152 * t160;
	t165 = t157 * t156;
	t164 = t157 * t159;
	t163 = t160 * t156;
	t162 = t160 * t159;
	t148 = cos(pkin(12)) * pkin(3) + pkin(2);
	t154 = qJ(3) + pkin(9);
	t161 = t148 * t156 - t154 * t159;
	t158 = cos(qJ(5));
	t155 = sin(qJ(5));
	t153 = cos(pkin(6));
	t151 = pkin(12) + qJ(4);
	t150 = cos(t151);
	t149 = sin(t151);
	t147 = sin(pkin(12)) * pkin(3) + pkin(8);
	t146 = t153 * t164 + t163;
	t145 = t153 * t163 + t164;
	t144 = -t153 * t162 + t165;
	t143 = t153 * t165 - t162;
	t142 = t148 * t159 + t154 * t156 + pkin(1);
	t141 = t153 * t149 + t150 * t169;
	t140 = t149 * t169 - t153 * t150;
	t139 = t152 * t147 - t161 * t153;
	t138 = -t143 * t150 + t149 * t168;
	t137 = t145 * t150 - t149 * t166;
	t136 = t145 * t149 + t150 * t166;
	t135 = t143 * t149 + t150 * t168;
	t1 = [t138 * t158 + t146 * t155, -t138 * t155 + t146 * t158, -t135, t138 * pkin(4) - t135 * pkin(10) + t139 * t157 + t142 * t160 + 0; t137 * t158 + t144 * t155, -t137 * t155 + t144 * t158, t136, t137 * pkin(4) + t136 * pkin(10) - t139 * t160 + t142 * t157 + 0; t141 * t158 - t155 * t167, -t141 * t155 - t158 * t167, t140, t141 * pkin(4) + t140 * pkin(10) + t147 * t153 + t161 * t152 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:03
	% EndTime: 2020-11-04 22:20:03
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (123->45), mult. (154->61), div. (0->0), fcn. (200->14), ass. (0->40)
	t209 = pkin(5) * sin(qJ(5));
	t183 = cos(pkin(12)) * pkin(3) + pkin(2);
	t195 = sin(qJ(2));
	t208 = t183 * t195;
	t191 = sin(pkin(6));
	t207 = t191 * t195;
	t196 = sin(qJ(1));
	t206 = t191 * t196;
	t197 = cos(qJ(2));
	t205 = t191 * t197;
	t198 = cos(qJ(1));
	t204 = t191 * t198;
	t203 = t196 * t195;
	t202 = t196 * t197;
	t201 = t198 * t195;
	t200 = t198 * t197;
	t199 = -pkin(11) - pkin(10);
	t193 = qJ(3) + pkin(9);
	t192 = cos(pkin(6));
	t190 = qJ(5) + qJ(6);
	t189 = pkin(12) + qJ(4);
	t188 = cos(t190);
	t187 = sin(t190);
	t186 = cos(t189);
	t185 = sin(t189);
	t184 = cos(qJ(5)) * pkin(5) + pkin(4);
	t182 = sin(pkin(12)) * pkin(3) + pkin(8);
	t181 = t192 * t202 + t201;
	t180 = t192 * t201 + t202;
	t179 = -t192 * t200 + t203;
	t178 = t192 * t203 - t200;
	t177 = t183 * t197 + t193 * t195 + pkin(1);
	t176 = t192 * t185 + t186 * t207;
	t175 = t185 * t207 - t192 * t186;
	t174 = t191 * t182 + (t193 * t197 - t208) * t192;
	t173 = -t178 * t186 + t185 * t206;
	t172 = t180 * t186 - t185 * t204;
	t171 = t180 * t185 + t186 * t204;
	t170 = t178 * t185 + t186 * t206;
	t1 = [t173 * t188 + t181 * t187, -t173 * t187 + t181 * t188, -t170, t170 * t199 + t173 * t184 + t174 * t196 + t177 * t198 + t181 * t209 + 0; t172 * t188 + t179 * t187, -t172 * t187 + t179 * t188, t171, -t171 * t199 + t172 * t184 - t174 * t198 + t177 * t196 + t179 * t209 + 0; t176 * t188 - t187 * t205, -t176 * t187 - t188 * t205, t175, -t175 * t199 + t176 * t184 + t182 * t192 + pkin(7) + 0 + (t208 + (-t193 - t209) * t197) * t191; 0, 0, 0, 1;];
	Tc_mdh = t1;
end