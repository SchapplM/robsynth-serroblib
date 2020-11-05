% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRP4 (for one body)
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

function Tc_mdh = S6PRPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:33
	% EndTime: 2020-11-04 21:02:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:33
	% EndTime: 2020-11-04 21:02:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t87 = cos(pkin(10));
	t86 = sin(pkin(10));
	t1 = [t87, -t86, 0, 0; t86, t87, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:33
	% EndTime: 2020-11-04 21:02:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t88 = sin(pkin(10));
	t89 = sin(pkin(6));
	t97 = t88 * t89;
	t90 = cos(pkin(10));
	t96 = t90 * t89;
	t91 = cos(pkin(6));
	t92 = sin(qJ(2));
	t95 = t91 * t92;
	t93 = cos(qJ(2));
	t94 = t91 * t93;
	t1 = [-t88 * t95 + t90 * t93, -t88 * t94 - t90 * t92, t97, t90 * pkin(1) + pkin(7) * t97 + 0; t88 * t93 + t90 * t95, -t88 * t92 + t90 * t94, -t96, t88 * pkin(1) - pkin(7) * t96 + 0; t89 * t92, t89 * t93, t91, t91 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:33
	% EndTime: 2020-11-04 21:02:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t101 = sin(pkin(10));
	t102 = sin(pkin(6));
	t114 = t101 * t102;
	t105 = cos(pkin(6));
	t113 = t101 * t105;
	t104 = cos(pkin(10));
	t112 = t102 * t104;
	t106 = sin(qJ(2));
	t111 = t102 * t106;
	t110 = t104 * t105;
	t109 = t105 * t106;
	t107 = cos(qJ(2));
	t108 = t105 * t107;
	t103 = cos(pkin(11));
	t100 = sin(pkin(11));
	t99 = t101 * t107 + t104 * t109;
	t98 = t101 * t109 - t104 * t107;
	t1 = [t100 * t114 - t98 * t103, t98 * t100 + t103 * t114, t101 * t108 + t104 * t106, (t104 * pkin(2) + qJ(3) * t113) * t107 + (-pkin(2) * t113 + t104 * qJ(3)) * t106 + pkin(7) * t114 + t104 * pkin(1) + 0; -t100 * t112 + t99 * t103, -t99 * t100 - t103 * t112, t101 * t106 - t104 * t108, (t101 * pkin(2) - qJ(3) * t110) * t107 + (pkin(2) * t110 + t101 * qJ(3)) * t106 - pkin(7) * t112 + t101 * pkin(1) + 0; t105 * t100 + t103 * t111, -t100 * t111 + t105 * t103, -t102 * t107, t105 * pkin(7) + qJ(1) + 0 + (pkin(2) * t106 - qJ(3) * t107) * t102; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:33
	% EndTime: 2020-11-04 21:02:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->31), mult. (72->49), div. (0->0), fcn. (93->10), ass. (0->23)
	t118 = cos(pkin(11)) * pkin(3) + pkin(2);
	t122 = sin(pkin(10));
	t136 = t122 * t118;
	t123 = sin(pkin(6));
	t135 = t122 * t123;
	t124 = cos(pkin(10));
	t134 = t123 * t124;
	t127 = sin(qJ(2));
	t133 = t123 * t127;
	t132 = t124 * t118;
	t125 = cos(pkin(6));
	t126 = qJ(3) + pkin(8);
	t131 = t125 * t126;
	t130 = t125 * t127;
	t128 = cos(qJ(2));
	t129 = t125 * t128;
	t121 = pkin(11) + qJ(4);
	t120 = cos(t121);
	t119 = sin(t121);
	t117 = sin(pkin(11)) * pkin(3) + pkin(7);
	t116 = t122 * t128 + t124 * t130;
	t115 = t122 * t130 - t124 * t128;
	t1 = [-t115 * t120 + t119 * t135, t115 * t119 + t120 * t135, t122 * t129 + t124 * t127, (t122 * t131 + t132) * t128 + (t124 * t126 - t125 * t136) * t127 + t117 * t135 + t124 * pkin(1) + 0; t116 * t120 - t119 * t134, -t116 * t119 - t120 * t134, t122 * t127 - t124 * t129, (-t124 * t131 + t136) * t128 + (t122 * t126 + t125 * t132) * t127 - t117 * t134 + t122 * pkin(1) + 0; t125 * t119 + t120 * t133, -t119 * t133 + t125 * t120, -t123 * t128, t117 * t125 + qJ(1) + 0 + (t118 * t127 - t126 * t128) * t123; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:33
	% EndTime: 2020-11-04 21:02:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (100->43), mult. (142->67), div. (0->0), fcn. (186->12), ass. (0->34)
	t148 = cos(pkin(11)) * pkin(3) + pkin(2);
	t152 = sin(pkin(10));
	t169 = t152 * t148;
	t153 = sin(pkin(6));
	t168 = t152 * t153;
	t154 = cos(pkin(10));
	t167 = t153 * t154;
	t158 = sin(qJ(2));
	t166 = t153 * t158;
	t160 = cos(qJ(2));
	t165 = t153 * t160;
	t164 = t154 * t148;
	t155 = cos(pkin(6));
	t156 = qJ(3) + pkin(8);
	t163 = t155 * t156;
	t162 = t155 * t158;
	t161 = t155 * t160;
	t159 = cos(qJ(5));
	t157 = sin(qJ(5));
	t151 = pkin(11) + qJ(4);
	t150 = cos(t151);
	t149 = sin(t151);
	t147 = sin(pkin(11)) * pkin(3) + pkin(7);
	t146 = t152 * t161 + t154 * t158;
	t145 = t152 * t160 + t154 * t162;
	t144 = t152 * t158 - t154 * t161;
	t143 = t152 * t162 - t154 * t160;
	t142 = t155 * t149 + t150 * t166;
	t141 = t149 * t166 - t155 * t150;
	t140 = -t143 * t150 + t149 * t168;
	t139 = t145 * t150 - t149 * t167;
	t138 = t145 * t149 + t150 * t167;
	t137 = t143 * t149 + t150 * t168;
	t1 = [t140 * t159 + t146 * t157, -t140 * t157 + t146 * t159, -t137, t140 * pkin(4) - t137 * pkin(9) + (t152 * t163 + t164) * t160 + (t154 * t156 - t155 * t169) * t158 + t147 * t168 + t154 * pkin(1) + 0; t139 * t159 + t144 * t157, -t139 * t157 + t144 * t159, t138, t139 * pkin(4) + t138 * pkin(9) + (-t154 * t163 + t169) * t160 + (t152 * t156 + t155 * t164) * t158 - t147 * t167 + t152 * pkin(1) + 0; t142 * t159 - t157 * t165, -t142 * t157 - t159 * t165, t141, t142 * pkin(4) + t141 * pkin(9) + t147 * t155 + qJ(1) + 0 + (t148 * t158 - t156 * t160) * t153; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:02:33
	% EndTime: 2020-11-04 21:02:33
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (138->49), mult. (204->73), div. (0->0), fcn. (268->12), ass. (0->40)
	t187 = cos(pkin(11)) * pkin(3) + pkin(2);
	t191 = sin(pkin(10));
	t208 = t191 * t187;
	t192 = sin(pkin(6));
	t207 = t191 * t192;
	t193 = cos(pkin(10));
	t206 = t192 * t193;
	t197 = sin(qJ(2));
	t205 = t192 * t197;
	t199 = cos(qJ(2));
	t204 = t192 * t199;
	t203 = t193 * t187;
	t194 = cos(pkin(6));
	t195 = qJ(3) + pkin(8);
	t202 = t194 * t195;
	t201 = t194 * t197;
	t200 = t194 * t199;
	t198 = cos(qJ(5));
	t196 = sin(qJ(5));
	t190 = pkin(11) + qJ(4);
	t189 = cos(t190);
	t188 = sin(t190);
	t186 = sin(pkin(11)) * pkin(3) + pkin(7);
	t185 = t191 * t200 + t193 * t197;
	t184 = t191 * t199 + t193 * t201;
	t183 = t191 * t197 - t193 * t200;
	t182 = t191 * t201 - t193 * t199;
	t181 = t194 * t188 + t189 * t205;
	t180 = t188 * t205 - t194 * t189;
	t179 = t181 * t198 - t196 * t204;
	t178 = t181 * t196 + t198 * t204;
	t177 = -t182 * t189 + t188 * t207;
	t176 = t184 * t189 - t188 * t206;
	t175 = t184 * t188 + t189 * t206;
	t174 = t182 * t188 + t189 * t207;
	t173 = t177 * t198 + t185 * t196;
	t172 = t177 * t196 - t185 * t198;
	t171 = t176 * t198 + t183 * t196;
	t170 = t176 * t196 - t183 * t198;
	t1 = [t173, -t174, t172, t173 * pkin(5) + t172 * qJ(6) + t177 * pkin(4) - t174 * pkin(9) + (t191 * t202 + t203) * t199 + (t193 * t195 - t194 * t208) * t197 + t186 * t207 + t193 * pkin(1) + 0; t171, t175, t170, t171 * pkin(5) + t170 * qJ(6) + t176 * pkin(4) + t175 * pkin(9) + (-t193 * t202 + t208) * t199 + (t191 * t195 + t194 * t203) * t197 - t186 * t206 + t191 * pkin(1) + 0; t179, t180, t178, t181 * pkin(4) + t179 * pkin(5) + t180 * pkin(9) + t178 * qJ(6) + t186 * t194 + qJ(1) + 0 + (t187 * t197 - t195 * t199) * t192; 0, 0, 0, 1;];
	Tc_mdh = t1;
end