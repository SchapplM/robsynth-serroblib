% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP9 (for one body)
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
% Datum: 2020-11-04 22:15
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:15:43
	% EndTime: 2020-11-04 22:15:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:15:43
	% EndTime: 2020-11-04 22:15:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t79 = cos(qJ(1));
	t78 = sin(qJ(1));
	t1 = [t79, -t78, 0, 0; t78, t79, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:15:43
	% EndTime: 2020-11-04 22:15:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t80 = sin(pkin(6));
	t83 = sin(qJ(1));
	t91 = t83 * t80;
	t82 = sin(qJ(2));
	t90 = t83 * t82;
	t84 = cos(qJ(2));
	t89 = t83 * t84;
	t85 = cos(qJ(1));
	t88 = t85 * t80;
	t87 = t85 * t82;
	t86 = t85 * t84;
	t81 = cos(pkin(6));
	t1 = [-t81 * t90 + t86, -t81 * t89 - t87, t91, t85 * pkin(1) + pkin(8) * t91 + 0; t81 * t87 + t89, t81 * t86 - t90, -t88, t83 * pkin(1) - pkin(8) * t88 + 0; t80 * t82, t80 * t84, t81, t81 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:15:43
	% EndTime: 2020-11-04 22:15:43
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->22), mult. (60->36), div. (0->0), fcn. (81->8), ass. (0->21)
	t100 = sin(qJ(2));
	t97 = sin(pkin(6));
	t111 = t100 * t97;
	t101 = sin(qJ(1));
	t110 = t101 * t97;
	t103 = cos(qJ(1));
	t109 = t103 * t97;
	t108 = t101 * t100;
	t102 = cos(qJ(2));
	t107 = t101 * t102;
	t106 = t103 * t100;
	t105 = t103 * t102;
	t104 = pkin(2) * t100 - qJ(3) * t102;
	t99 = cos(pkin(6));
	t98 = cos(pkin(11));
	t96 = sin(pkin(11));
	t95 = t102 * pkin(2) + t100 * qJ(3) + pkin(1);
	t94 = -t99 * t108 + t105;
	t93 = t99 * t106 + t107;
	t92 = t97 * pkin(8) - t104 * t99;
	t1 = [t96 * t110 + t94 * t98, t98 * t110 - t94 * t96, t99 * t107 + t106, t92 * t101 + t95 * t103 + 0; -t96 * t109 + t93 * t98, -t98 * t109 - t93 * t96, -t99 * t105 + t108, t95 * t101 - t92 * t103 + 0; t98 * t111 + t99 * t96, -t96 * t111 + t99 * t98, -t97 * t102, t99 * pkin(8) + t104 * t97 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:15:43
	% EndTime: 2020-11-04 22:15:43
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (54->26), mult. (68->38), div. (0->0), fcn. (89->10), ass. (0->25)
	t121 = sin(pkin(6));
	t124 = sin(qJ(2));
	t135 = t121 * t124;
	t125 = sin(qJ(1));
	t134 = t121 * t125;
	t127 = cos(qJ(1));
	t133 = t121 * t127;
	t132 = t125 * t124;
	t126 = cos(qJ(2));
	t131 = t125 * t126;
	t130 = t127 * t124;
	t129 = t127 * t126;
	t117 = cos(pkin(11)) * pkin(3) + pkin(2);
	t123 = qJ(3) + pkin(9);
	t128 = t117 * t124 - t123 * t126;
	t122 = cos(pkin(6));
	t120 = pkin(11) + qJ(4);
	t119 = cos(t120);
	t118 = sin(t120);
	t116 = sin(pkin(11)) * pkin(3) + pkin(8);
	t115 = -t122 * t132 + t129;
	t114 = t122 * t130 + t131;
	t113 = t117 * t126 + t123 * t124 + pkin(1);
	t112 = t121 * t116 - t128 * t122;
	t1 = [t115 * t119 + t118 * t134, -t115 * t118 + t119 * t134, t122 * t131 + t130, t112 * t125 + t113 * t127 + 0; t114 * t119 - t118 * t133, -t114 * t118 - t119 * t133, -t122 * t129 + t132, -t112 * t127 + t113 * t125 + 0; t122 * t118 + t119 * t135, -t118 * t135 + t122 * t119, -t121 * t126, t116 * t122 + t128 * t121 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:15:43
	% EndTime: 2020-11-04 22:15:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (100->38), mult. (138->56), div. (0->0), fcn. (182->12), ass. (0->36)
	t153 = sin(pkin(6));
	t157 = sin(qJ(2));
	t170 = t153 * t157;
	t158 = sin(qJ(1));
	t169 = t153 * t158;
	t160 = cos(qJ(2));
	t168 = t153 * t160;
	t161 = cos(qJ(1));
	t167 = t153 * t161;
	t166 = t158 * t157;
	t165 = t158 * t160;
	t164 = t161 * t157;
	t163 = t161 * t160;
	t149 = cos(pkin(11)) * pkin(3) + pkin(2);
	t155 = qJ(3) + pkin(9);
	t162 = t149 * t157 - t155 * t160;
	t159 = cos(qJ(5));
	t156 = sin(qJ(5));
	t154 = cos(pkin(6));
	t152 = pkin(11) + qJ(4);
	t151 = cos(t152);
	t150 = sin(t152);
	t148 = sin(pkin(11)) * pkin(3) + pkin(8);
	t147 = -t154 * t166 + t163;
	t146 = t154 * t165 + t164;
	t145 = t154 * t164 + t165;
	t144 = -t154 * t163 + t166;
	t143 = t149 * t160 + t155 * t157 + pkin(1);
	t142 = t154 * t150 + t151 * t170;
	t141 = t150 * t170 - t154 * t151;
	t140 = t153 * t148 - t162 * t154;
	t139 = -t147 * t150 + t151 * t169;
	t138 = t147 * t151 + t150 * t169;
	t137 = t145 * t151 - t150 * t167;
	t136 = t145 * t150 + t151 * t167;
	t1 = [t138 * t159 + t146 * t156, -t138 * t156 + t146 * t159, -t139, t138 * pkin(4) - t139 * pkin(10) + t140 * t158 + t143 * t161 + 0; t137 * t159 + t144 * t156, -t137 * t156 + t144 * t159, t136, t137 * pkin(4) + t136 * pkin(10) - t140 * t161 + t143 * t158 + 0; t142 * t159 - t156 * t168, -t142 * t156 - t159 * t168, t141, t142 * pkin(4) + t141 * pkin(10) + t148 * t154 + t162 * t153 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:15:43
	% EndTime: 2020-11-04 22:15:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (111->44), mult. (154->61), div. (0->0), fcn. (200->12), ass. (0->40)
	t190 = cos(pkin(6));
	t197 = cos(qJ(2));
	t198 = cos(qJ(1));
	t199 = t198 * t197;
	t194 = sin(qJ(2));
	t195 = sin(qJ(1));
	t202 = t195 * t194;
	t179 = -t190 * t199 + t202;
	t193 = sin(qJ(5));
	t209 = t179 * t193;
	t200 = t198 * t194;
	t201 = t195 * t197;
	t181 = t190 * t201 + t200;
	t208 = t181 * t193;
	t184 = cos(pkin(11)) * pkin(3) + pkin(2);
	t207 = t184 * t194;
	t189 = sin(pkin(6));
	t206 = t189 * t194;
	t205 = t189 * t195;
	t204 = t189 * t197;
	t203 = t189 * t198;
	t196 = cos(qJ(5));
	t192 = qJ(3) + pkin(9);
	t191 = -qJ(6) - pkin(10);
	t188 = pkin(11) + qJ(4);
	t187 = cos(t188);
	t186 = sin(t188);
	t185 = t196 * pkin(5) + pkin(4);
	t183 = sin(pkin(11)) * pkin(3) + pkin(8);
	t182 = -t190 * t202 + t199;
	t180 = t190 * t200 + t201;
	t178 = t184 * t197 + t192 * t194 + pkin(1);
	t177 = t190 * t186 + t187 * t206;
	t176 = t186 * t206 - t190 * t187;
	t175 = t189 * t183 + (t192 * t197 - t207) * t190;
	t174 = -t182 * t186 + t187 * t205;
	t173 = t182 * t187 + t186 * t205;
	t172 = t180 * t187 - t186 * t203;
	t171 = t180 * t186 + t187 * t203;
	t1 = [t173 * t196 + t208, -t173 * t193 + t181 * t196, -t174, pkin(5) * t208 + t173 * t185 + t174 * t191 + t175 * t195 + t178 * t198 + 0; t172 * t196 + t209, -t172 * t193 + t179 * t196, t171, pkin(5) * t209 - t171 * t191 + t172 * t185 - t175 * t198 + t178 * t195 + 0; t177 * t196 - t193 * t204, -t177 * t193 - t196 * t204, t176, -t176 * t191 + t177 * t185 + t183 * t190 + pkin(7) + 0 + (t207 + (-pkin(5) * t193 - t192) * t197) * t189; 0, 0, 0, 1;];
	Tc_mdh = t1;
end