% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:02
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:09
	% EndTime: 2020-11-04 21:02:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:09
	% EndTime: 2020-11-04 21:02:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t84 = cos(pkin(10));
	t83 = sin(pkin(10));
	t1 = [t84, -t83, 0, 0; t83, t84, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:09
	% EndTime: 2020-11-04 21:02:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t85 = sin(pkin(10));
	t86 = sin(pkin(6));
	t94 = t85 * t86;
	t87 = cos(pkin(10));
	t93 = t87 * t86;
	t88 = cos(pkin(6));
	t89 = sin(qJ(2));
	t92 = t88 * t89;
	t90 = cos(qJ(2));
	t91 = t88 * t90;
	t1 = [-t85 * t92 + t87 * t90, -t85 * t91 - t87 * t89, t94, t87 * pkin(1) + pkin(7) * t94 + 0; t85 * t90 + t87 * t92, -t85 * t89 + t87 * t91, -t93, t85 * pkin(1) - pkin(7) * t93 + 0; t86 * t89, t86 * t90, t88, t88 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:09
	% EndTime: 2020-11-04 21:02:09
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (29->27), mult. (64->47), div. (0->0), fcn. (85->8), ass. (0->19)
	t98 = sin(pkin(10));
	t112 = t98 * pkin(2);
	t99 = sin(pkin(6));
	t111 = t98 * t99;
	t101 = cos(pkin(10));
	t110 = t101 * t99;
	t103 = sin(qJ(2));
	t109 = t103 * t99;
	t108 = t98 * qJ(3);
	t102 = cos(pkin(6));
	t107 = t101 * t102;
	t106 = t102 * t103;
	t104 = cos(qJ(2));
	t105 = t102 * t104;
	t100 = cos(pkin(11));
	t97 = sin(pkin(11));
	t96 = t101 * t106 + t98 * t104;
	t95 = -t101 * t104 + t98 * t106;
	t1 = [-t95 * t100 + t97 * t111, t100 * t111 + t95 * t97, t101 * t103 + t98 * t105, (t101 * pkin(2) + t102 * t108) * t104 + (t101 * qJ(3) - t102 * t112) * t103 + pkin(7) * t111 + t101 * pkin(1) + 0; t96 * t100 - t97 * t110, -t100 * t110 - t96 * t97, -t101 * t105 + t98 * t103, (-qJ(3) * t107 + t112) * t104 + (pkin(2) * t107 + t108) * t103 - pkin(7) * t110 + t98 * pkin(1) + 0; t100 * t109 + t102 * t97, t102 * t100 - t97 * t109, -t99 * t104, t102 * pkin(7) + qJ(1) + 0 + (pkin(2) * t103 - qJ(3) * t104) * t99; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:09
	% EndTime: 2020-11-04 21:02:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->31), mult. (72->49), div. (0->0), fcn. (93->10), ass. (0->23)
	t116 = cos(pkin(11)) * pkin(3) + pkin(2);
	t120 = sin(pkin(10));
	t134 = t120 * t116;
	t121 = sin(pkin(6));
	t133 = t120 * t121;
	t122 = cos(pkin(10));
	t132 = t121 * t122;
	t125 = sin(qJ(2));
	t131 = t121 * t125;
	t130 = t122 * t116;
	t123 = cos(pkin(6));
	t124 = qJ(3) + pkin(8);
	t129 = t123 * t124;
	t128 = t123 * t125;
	t126 = cos(qJ(2));
	t127 = t123 * t126;
	t119 = pkin(11) + qJ(4);
	t118 = cos(t119);
	t117 = sin(t119);
	t115 = sin(pkin(11)) * pkin(3) + pkin(7);
	t114 = t120 * t126 + t122 * t128;
	t113 = t120 * t128 - t122 * t126;
	t1 = [-t113 * t118 + t117 * t133, t113 * t117 + t118 * t133, t120 * t127 + t122 * t125, (t120 * t129 + t130) * t126 + (t122 * t124 - t123 * t134) * t125 + t115 * t133 + t122 * pkin(1) + 0; t114 * t118 - t117 * t132, -t114 * t117 - t118 * t132, t120 * t125 - t122 * t127, (-t122 * t129 + t134) * t126 + (t120 * t124 + t123 * t130) * t125 - t115 * t132 + t120 * pkin(1) + 0; t123 * t117 + t118 * t131, -t117 * t131 + t123 * t118, -t121 * t126, t115 * t123 + qJ(1) + 0 + (t116 * t125 - t124 * t126) * t121; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:09
	% EndTime: 2020-11-04 21:02:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (100->43), mult. (142->67), div. (0->0), fcn. (186->12), ass. (0->34)
	t146 = cos(pkin(11)) * pkin(3) + pkin(2);
	t150 = sin(pkin(10));
	t167 = t150 * t146;
	t151 = sin(pkin(6));
	t166 = t150 * t151;
	t152 = cos(pkin(10));
	t165 = t151 * t152;
	t156 = sin(qJ(2));
	t164 = t151 * t156;
	t158 = cos(qJ(2));
	t163 = t151 * t158;
	t162 = t152 * t146;
	t153 = cos(pkin(6));
	t154 = qJ(3) + pkin(8);
	t161 = t153 * t154;
	t160 = t153 * t156;
	t159 = t153 * t158;
	t157 = cos(qJ(5));
	t155 = sin(qJ(5));
	t149 = pkin(11) + qJ(4);
	t148 = cos(t149);
	t147 = sin(t149);
	t145 = sin(pkin(11)) * pkin(3) + pkin(7);
	t144 = t150 * t159 + t152 * t156;
	t143 = t150 * t158 + t152 * t160;
	t142 = t150 * t156 - t152 * t159;
	t141 = t150 * t160 - t152 * t158;
	t140 = t153 * t147 + t148 * t164;
	t139 = t147 * t164 - t153 * t148;
	t138 = -t141 * t148 + t147 * t166;
	t137 = t143 * t148 - t147 * t165;
	t136 = t143 * t147 + t148 * t165;
	t135 = t141 * t147 + t148 * t166;
	t1 = [t138 * t157 + t144 * t155, -t138 * t155 + t144 * t157, -t135, t138 * pkin(4) - t135 * pkin(9) + (t150 * t161 + t162) * t158 + (t152 * t154 - t153 * t167) * t156 + t145 * t166 + t152 * pkin(1) + 0; t137 * t157 + t142 * t155, -t137 * t155 + t142 * t157, t136, t137 * pkin(4) + t136 * pkin(9) + (-t152 * t161 + t167) * t158 + (t150 * t154 + t153 * t162) * t156 - t145 * t165 + t150 * pkin(1) + 0; t140 * t157 - t155 * t163, -t140 * t155 - t157 * t163, t139, t140 * pkin(4) + t139 * pkin(9) + t145 * t153 + qJ(1) + 0 + (t146 * t156 - t154 * t158) * t151; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:09
	% EndTime: 2020-11-04 21:02:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (111->48), mult. (158->71), div. (0->0), fcn. (204->12), ass. (0->38)
	t184 = sin(pkin(10));
	t186 = cos(pkin(10));
	t191 = sin(qJ(2));
	t187 = cos(pkin(6));
	t193 = cos(qJ(2));
	t194 = t187 * t193;
	t175 = t184 * t191 - t186 * t194;
	t190 = sin(qJ(5));
	t204 = t175 * t190;
	t177 = t184 * t194 + t186 * t191;
	t203 = t177 * t190;
	t179 = cos(pkin(11)) * pkin(3) + pkin(2);
	t202 = t184 * t179;
	t185 = sin(pkin(6));
	t201 = t184 * t185;
	t200 = t185 * t186;
	t199 = t185 * t191;
	t198 = t185 * t193;
	t197 = t186 * t179;
	t189 = qJ(3) + pkin(8);
	t196 = t187 * t189;
	t195 = t187 * t191;
	t192 = cos(qJ(5));
	t188 = -qJ(6) - pkin(9);
	t183 = pkin(11) + qJ(4);
	t182 = cos(t183);
	t181 = sin(t183);
	t180 = t192 * pkin(5) + pkin(4);
	t178 = sin(pkin(11)) * pkin(3) + pkin(7);
	t176 = t184 * t193 + t186 * t195;
	t174 = t184 * t195 - t186 * t193;
	t173 = t187 * t181 + t182 * t199;
	t172 = t181 * t199 - t187 * t182;
	t171 = -t174 * t182 + t181 * t201;
	t170 = t176 * t182 - t181 * t200;
	t169 = t176 * t181 + t182 * t200;
	t168 = t174 * t181 + t182 * t201;
	t1 = [t171 * t192 + t203, -t171 * t190 + t177 * t192, -t168, t171 * t180 + t168 * t188 + pkin(5) * t203 + (t184 * t196 + t197) * t193 + (t186 * t189 - t187 * t202) * t191 + t178 * t201 + t186 * pkin(1) + 0; t170 * t192 + t204, -t170 * t190 + t175 * t192, t169, t170 * t180 - t169 * t188 + pkin(5) * t204 + (-t186 * t196 + t202) * t193 + (t184 * t189 + t187 * t197) * t191 - t178 * t200 + t184 * pkin(1) + 0; t173 * t192 - t190 * t198, -t173 * t190 - t192 * t198, t172, -t172 * t188 + t173 * t180 + t178 * t187 + qJ(1) + 0 + (t179 * t191 + (-pkin(5) * t190 - t189) * t193) * t185; 0, 0, 0, 1;];
	Tc_mdh = t1;
end