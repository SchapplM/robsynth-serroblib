% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:16
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:06
	% EndTime: 2020-11-04 22:16:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:06
	% EndTime: 2020-11-04 22:16:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t82 = cos(qJ(1));
	t81 = sin(qJ(1));
	t1 = [t82, -t81, 0, 0; t81, t82, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:06
	% EndTime: 2020-11-04 22:16:06
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
	% StartTime: 2020-11-04 22:16:06
	% EndTime: 2020-11-04 22:16:06
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
	% StartTime: 2020-11-04 22:16:06
	% EndTime: 2020-11-04 22:16:06
	% DurationCPUTime: 0.10s
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
	t115 = t124 * t119 - t125 * t131;
	t1 = [t118 * t122 + t121 * t137, -t118 * t121 + t122 * t137, t125 * t134 + t133, t115 * t128 + t116 * t130 + 0; t117 * t122 - t121 * t136, -t117 * t121 - t122 * t136, -t125 * t132 + t135, -t115 * t130 + t116 * t128 + 0; t125 * t121 + t122 * t138, -t121 * t138 + t125 * t122, -t124 * t129, t119 * t125 + t124 * t131 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:06
	% EndTime: 2020-11-04 22:16:06
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (100->38), mult. (138->56), div. (0->0), fcn. (182->12), ass. (0->36)
	t156 = sin(pkin(6));
	t160 = sin(qJ(2));
	t173 = t156 * t160;
	t161 = sin(qJ(1));
	t172 = t156 * t161;
	t163 = cos(qJ(2));
	t171 = t156 * t163;
	t164 = cos(qJ(1));
	t170 = t156 * t164;
	t169 = t161 * t160;
	t168 = t161 * t163;
	t167 = t164 * t160;
	t166 = t164 * t163;
	t152 = cos(pkin(11)) * pkin(3) + pkin(2);
	t158 = qJ(3) + pkin(9);
	t165 = t152 * t160 - t158 * t163;
	t162 = cos(qJ(5));
	t159 = sin(qJ(5));
	t157 = cos(pkin(6));
	t155 = pkin(11) + qJ(4);
	t154 = cos(t155);
	t153 = sin(t155);
	t151 = sin(pkin(11)) * pkin(3) + pkin(8);
	t150 = -t157 * t169 + t166;
	t149 = t157 * t168 + t167;
	t148 = t157 * t167 + t168;
	t147 = -t157 * t166 + t169;
	t146 = t152 * t163 + t158 * t160 + pkin(1);
	t145 = t157 * t153 + t154 * t173;
	t144 = t153 * t173 - t157 * t154;
	t143 = t156 * t151 - t165 * t157;
	t142 = -t150 * t153 + t154 * t172;
	t141 = t150 * t154 + t153 * t172;
	t140 = t148 * t154 - t153 * t170;
	t139 = t148 * t153 + t154 * t170;
	t1 = [t141 * t162 + t149 * t159, -t141 * t159 + t149 * t162, -t142, t141 * pkin(4) - t142 * pkin(10) + t143 * t161 + t146 * t164 + 0; t140 * t162 + t147 * t159, -t140 * t159 + t147 * t162, t139, t140 * pkin(4) + t139 * pkin(10) - t143 * t164 + t146 * t161 + 0; t145 * t162 - t159 * t171, -t145 * t159 - t162 * t171, t144, t145 * pkin(4) + t144 * pkin(10) + t151 * t157 + t165 * t156 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:06
	% EndTime: 2020-11-04 22:16:06
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (138->44), mult. (200->62), div. (0->0), fcn. (264->12), ass. (0->42)
	t197 = sin(pkin(6));
	t201 = sin(qJ(2));
	t214 = t197 * t201;
	t202 = sin(qJ(1));
	t213 = t197 * t202;
	t204 = cos(qJ(2));
	t212 = t197 * t204;
	t205 = cos(qJ(1));
	t211 = t197 * t205;
	t210 = t202 * t201;
	t209 = t202 * t204;
	t208 = t205 * t201;
	t207 = t205 * t204;
	t193 = cos(pkin(11)) * pkin(3) + pkin(2);
	t199 = qJ(3) + pkin(9);
	t206 = t193 * t201 - t199 * t204;
	t203 = cos(qJ(5));
	t200 = sin(qJ(5));
	t198 = cos(pkin(6));
	t196 = pkin(11) + qJ(4);
	t195 = cos(t196);
	t194 = sin(t196);
	t192 = sin(pkin(11)) * pkin(3) + pkin(8);
	t191 = -t198 * t210 + t207;
	t190 = t198 * t209 + t208;
	t189 = t198 * t208 + t209;
	t188 = -t198 * t207 + t210;
	t187 = t193 * t204 + t199 * t201 + pkin(1);
	t186 = t198 * t194 + t195 * t214;
	t185 = t194 * t214 - t198 * t195;
	t184 = t197 * t192 - t206 * t198;
	t183 = -t191 * t194 + t195 * t213;
	t182 = t191 * t195 + t194 * t213;
	t181 = t189 * t195 - t194 * t211;
	t180 = t189 * t194 + t195 * t211;
	t179 = t186 * t203 - t200 * t212;
	t178 = t186 * t200 + t203 * t212;
	t177 = t182 * t203 + t190 * t200;
	t176 = t182 * t200 - t190 * t203;
	t175 = t181 * t203 + t188 * t200;
	t174 = t181 * t200 - t188 * t203;
	t1 = [t177, -t183, t176, t182 * pkin(4) + t177 * pkin(5) - t183 * pkin(10) + t176 * qJ(6) + t184 * t202 + t187 * t205 + 0; t175, t180, t174, t181 * pkin(4) + t175 * pkin(5) + t180 * pkin(10) + t174 * qJ(6) - t184 * t205 + t187 * t202 + 0; t179, t185, t178, t186 * pkin(4) + t179 * pkin(5) + t185 * pkin(10) + t178 * qJ(6) + t192 * t198 + t206 * t197 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end