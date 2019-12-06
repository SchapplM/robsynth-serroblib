% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:47
% EndTime: 2019-12-05 16:33:00
% DurationCPUTime: 2.95s
% Computational Cost: add. (2959->326), mult. (7848->499), div. (0->0), fcn. (5875->10), ass. (0->179)
t129 = sin(pkin(10));
t134 = sin(qJ(3));
t188 = qJD(3) * t134;
t181 = pkin(7) * t188;
t116 = t129 * t181;
t131 = cos(pkin(10));
t135 = sin(qJ(2));
t130 = sin(pkin(5));
t193 = qJD(1) * t130;
t137 = cos(qJ(3));
t138 = cos(qJ(2));
t199 = t137 * t138;
t155 = pkin(3) * t134 - qJ(4) * t137;
t85 = qJD(3) * t155 - t134 * qJD(4);
t218 = t131 * t85 + t116 - (-t129 * t199 + t131 * t135) * t193;
t232 = t129 * t85 - (t129 * t135 + t131 * t199) * t193;
t133 = sin(qJ(5));
t136 = cos(qJ(5));
t184 = t131 * qJD(3);
t189 = qJD(2) * t134;
t97 = -t129 * t189 + t184;
t171 = t131 * t189;
t185 = t129 * qJD(3);
t98 = t171 + t185;
t47 = t133 * t98 - t136 * t97;
t231 = t47 ^ 2;
t50 = t133 * t97 + t136 * t98;
t230 = t50 ^ 2;
t202 = t131 * t137;
t150 = pkin(4) * t134 - pkin(8) * t202;
t147 = t150 * qJD(3);
t229 = -t147 - t218;
t203 = t131 * t134;
t207 = t129 * t137;
t228 = (-pkin(7) * t203 - pkin(8) * t207) * qJD(3) + t232;
t183 = t137 * qJD(2);
t121 = -qJD(5) + t183;
t227 = t47 * t121;
t186 = qJD(5) * t136;
t200 = t133 * t129;
t226 = -qJD(5) * t200 + t131 * t186;
t182 = qJD(2) * qJD(3);
t168 = t137 * t182;
t158 = t129 * t168;
t160 = t136 * t137 * t184;
t22 = t133 * (qJD(5) * t98 + t158) - qJD(2) * t160 - t97 * t186;
t172 = t138 * t193;
t159 = qJD(2) * t172;
t132 = cos(pkin(5));
t201 = t132 * t137;
t169 = qJD(1) * t201;
t196 = qJD(3) * t169 + t137 * t159;
t173 = t135 * t193;
t106 = qJD(2) * pkin(7) + t173;
t94 = t134 * t106;
t42 = (qJD(4) - t94) * qJD(3) + t196;
t55 = (t85 + t173) * qJD(2);
t14 = -t129 * t42 + t131 * t55;
t12 = qJD(2) * t147 + t14;
t15 = t129 * t55 + t131 * t42;
t13 = -pkin(8) * t158 + t15;
t192 = qJD(1) * t134;
t170 = t132 * t192;
t74 = t137 * t106 + t170;
t70 = qJD(3) * qJ(4) + t74;
t108 = -t137 * pkin(3) - t134 * qJ(4) - pkin(2);
t75 = qJD(2) * t108 - t172;
t25 = -t129 * t70 + t131 * t75;
t16 = -pkin(4) * t183 - t98 * pkin(8) + t25;
t26 = t129 * t75 + t131 * t70;
t17 = t97 * pkin(8) + t26;
t153 = t133 * t17 - t136 * t16;
t1 = -qJD(5) * t153 + t133 * t12 + t136 * t13;
t96 = t131 * t108;
t53 = -pkin(8) * t203 + t96 + (-pkin(7) * t129 - pkin(4)) * t137;
t77 = pkin(7) * t202 + t129 * t108;
t63 = -t129 * t134 * pkin(8) + t77;
t19 = t133 * t53 + t136 * t63;
t225 = qJD(5) * t19 + t228 * t133 + t229 * t136;
t18 = -t133 * t63 + t136 * t53;
t224 = -qJD(5) * t18 + t229 * t133 - t228 * t136;
t101 = t136 * t129 + t133 * t131;
t102 = t155 * qJD(2);
t73 = -t94 + t169;
t37 = t131 * t102 - t129 * t73;
t24 = qJD(2) * t150 + t37;
t178 = t129 * t183;
t38 = t129 * t102 + t131 * t73;
t29 = -pkin(8) * t178 + t38;
t220 = pkin(8) + qJ(4);
t111 = t220 * t129;
t112 = t220 * t131;
t65 = -t133 * t111 + t136 * t112;
t223 = -qJD(4) * t101 - qJD(5) * t65 + t133 * t29 - t136 * t24;
t190 = qJD(2) * t130;
t176 = t138 * t190;
t187 = qJD(3) * t137;
t45 = t106 * t187 + (qJD(3) * t132 + t176) * t192;
t206 = t130 * t135;
t89 = t134 * t206 - t201;
t222 = t45 * t89;
t221 = t50 * t47;
t100 = -t136 * t131 + t200;
t64 = -t136 * t111 - t133 * t112;
t219 = qJD(4) * t100 - qJD(5) * t64 + t133 * t24 + t136 * t29;
t164 = t131 * t181;
t217 = -t164 + t232;
t148 = t101 * t137;
t88 = t101 * qJD(5);
t216 = qJD(2) * t148 - t88;
t215 = -t100 * t183 - t226;
t214 = qJD(2) * pkin(2);
t213 = t131 * t97;
t212 = t45 * t129;
t211 = t45 * t131;
t210 = t45 * t134;
t128 = t137 ^ 2;
t140 = qJD(2) ^ 2;
t208 = t128 * t140;
t205 = t130 * t138;
t204 = t130 * t140;
t139 = qJD(3) ^ 2;
t198 = t139 * t134;
t197 = t139 * t137;
t127 = t134 ^ 2;
t195 = t127 - 0.2e1 * t128;
t194 = t127 - t128;
t191 = qJD(2) * t129;
t180 = t135 * t204;
t179 = pkin(4) * t129 + pkin(7);
t177 = t135 * t190;
t167 = t134 * t182;
t166 = -qJD(3) * pkin(3) + qJD(4);
t165 = -t97 + t184;
t163 = t134 * t172;
t162 = t134 * t176;
t161 = t137 * t176;
t157 = t137 * t167;
t107 = -t172 - t214;
t156 = -t107 - t172;
t8 = t133 * t16 + t136 * t17;
t90 = t132 * t134 + t137 * t206;
t59 = -t90 * t129 - t131 * t205;
t60 = -t129 * t205 + t90 * t131;
t20 = -t133 * t60 + t136 * t59;
t21 = t133 * t59 + t136 * t60;
t152 = qJD(2) * (-t98 + t185);
t151 = qJD(2) * t165;
t149 = t137 * t152;
t145 = qJD(3) * t148;
t66 = t166 - t73;
t144 = qJD(3) * (-t156 - t214);
t143 = -qJ(4) * t188 + (t166 - t66) * t137;
t2 = -qJD(5) * t8 + t136 * t12 - t133 * t13;
t44 = -t106 * t188 + t196;
t141 = t210 + t44 * t137 + (-t134 * t74 - t137 * t73) * qJD(3);
t30 = pkin(4) * t158 + t45;
t23 = qJD(2) * t145 + qJD(5) * t50;
t126 = t131 ^ 2;
t125 = t129 ^ 2;
t124 = -t131 * pkin(4) - pkin(3);
t119 = t134 * t140 * t137;
t114 = -0.2e1 * t157;
t103 = t179 * t134;
t93 = t179 * t187;
t83 = t100 * t134;
t82 = t101 * t134;
t76 = -pkin(7) * t207 + t96;
t62 = qJD(3) * t90 + t162;
t61 = -qJD(3) * t89 + t161;
t54 = t170 + (pkin(4) * t191 + t106) * t137;
t41 = t226 * t134 + t145;
t40 = t133 * t137 * t185 + t134 * t88 - t160;
t39 = -t97 * pkin(4) + t66;
t36 = t129 * t177 + t61 * t131;
t35 = -t61 * t129 + t131 * t177;
t4 = -qJD(5) * t21 - t133 * t36 + t136 * t35;
t3 = qJD(5) * t20 + t133 * t35 + t136 * t36;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, -t138 * t204, 0, 0, 0, 0, 0, 0, 0, 0, -t137 * t180 + (-t62 - t162) * qJD(3), t134 * t180 + (-t61 - t161) * qJD(3), (t134 * t62 + t137 * t61 + (-t134 * t90 + t137 * t89) * qJD(3)) * qJD(2), t44 * t90 + t222 + t74 * t61 - t73 * t62 + (t107 - t172) * t177, 0, 0, 0, 0, 0, 0, -t62 * t97 + (-t137 * t35 + (t134 * t59 + t89 * t207) * qJD(3)) * qJD(2), t62 * t98 + (t137 * t36 + (-t134 * t60 + t89 * t202) * qJD(3)) * qJD(2), -t35 * t98 + t36 * t97 + (-t129 * t60 - t131 * t59) * t168, t14 * t59 + t15 * t60 + t25 * t35 + t26 * t36 + t66 * t62 + t222, 0, 0, 0, 0, 0, 0, -t4 * t121 + t167 * t20 + t89 * t23 + t62 * t47, t3 * t121 - t167 * t21 - t89 * t22 + t62 * t50, t20 * t22 - t21 * t23 - t3 * t47 - t4 * t50, t1 * t21 - t153 * t4 + t2 * t20 + t8 * t3 + t30 * t89 + t39 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t157, -0.2e1 * t194 * t182, t197, t114, -t198, 0, -pkin(7) * t197 + t134 * t144, pkin(7) * t198 + t137 * t144, (-t127 - t128) * t159 + t141, ((t134 * t73 - t137 * t74) * t138 + (-t107 - t214) * t135) * t193 + t141 * pkin(7), (t126 * t189 + t131 * t98) * t187, (t213 + (-t98 - 0.2e1 * t171) * t129) * t187, (t195 * t131 * qJD(2) + t134 * t98) * qJD(3), (t125 * t189 - t129 * t97) * t187, (t134 * t97 - t195 * t191) * qJD(3), t114, (t97 * t172 + t212 + (qJD(2) * t76 + t25) * qJD(3)) * t134 + (-t14 + (-pkin(7) * t97 + t129 * t66) * qJD(3) + (t116 - t218) * qJD(2)) * t137, (-t98 * t172 + t211 + (-qJD(2) * t77 - t26) * qJD(3)) * t134 + (t15 + (pkin(7) * t98 + t131 * t66) * qJD(3) + (t164 + t217) * qJD(2)) * t137, -t218 * t98 + t217 * t97 + (-t129 * t15 - t131 * t14) * t134 + (-t129 * t26 - t131 * t25 + (-t129 * t77 - t131 * t76) * qJD(2)) * t187, -t66 * t163 + t14 * t76 + t15 * t77 + t217 * t26 + t218 * t25 + (t66 * t187 + t210) * pkin(7), t22 * t83 - t50 * t40, t22 * t82 + t83 * t23 + t40 * t47 - t50 * t41, t40 * t121 + t22 * t137 + (-qJD(2) * t83 + t50) * t188, t23 * t82 + t47 * t41, t41 * t121 + t23 * t137 + (-qJD(2) * t82 - t47) * t188, (-t121 - t183) * t188, t103 * t23 - t2 * t137 + t30 * t82 + t39 * t41 + t93 * t47 + t225 * t121 + (-t47 * t172 + (qJD(2) * t18 - t153) * qJD(3)) * t134, t1 * t137 - t103 * t22 - t30 * t83 - t39 * t40 + t93 * t50 - t224 * t121 + (-t50 * t172 + (-qJD(2) * t19 - t8) * qJD(3)) * t134, -t1 * t82 - t153 * t40 + t18 * t22 - t19 * t23 + t2 * t83 + t224 * t47 + t225 * t50 - t8 * t41, t1 * t19 + t30 * t103 + t2 * t18 - t224 * t8 + t225 * t153 + (t93 - t163) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t194 * t140, 0, t119, 0, 0, t156 * t189, -t107 * t183 + (t73 + t94) * qJD(3) - t196, 0, 0, t131 * t149, (t129 * t98 - t213 + (-t125 + t126) * qJD(3)) * t183, t131 * t208 + t134 * t152, -t165 * t178, -t129 * t208 + t134 * t151, t119, -t211 + t74 * t97 + (t129 * t143 - t134 * t25 + t137 * t37) * qJD(2), t212 - t74 * t98 + (t131 * t143 + t134 * t26 - t137 * t38) * qJD(2), t37 * t98 - t38 * t97 + (qJD(4) * t97 + t25 * t183 + t15) * t131 + (qJD(4) * t98 + t26 * t183 - t14) * t129, -t45 * pkin(3) - t25 * t37 - t26 * t38 - t66 * t74 + (-t129 * t25 + t131 * t26) * qJD(4) + (-t14 * t129 + t15 * t131) * qJ(4), -t22 * t101 - t215 * t50, t22 * t100 - t101 * t23 + t215 * t47 + t216 * t50, t215 * t121 + (qJD(3) * t101 - t50) * t189, t23 * t100 - t216 * t47, -t216 * t121 + (-qJD(3) * t100 + t47) * t189, t121 * t189, t30 * t100 + t124 * t23 - t54 * t47 - t216 * t39 - t223 * t121 + (qJD(3) * t64 + t153) * t189, t30 * t101 - t124 * t22 - t54 * t50 - t215 * t39 - t219 * t121 + (-qJD(3) * t65 + t8) * t189, -t1 * t100 - t2 * t101 - t153 * t215 + t216 * t8 + t219 * t47 + t64 * t22 - t223 * t50 - t65 * t23, t1 * t65 + t30 * t124 - t153 * t223 + t2 * t64 - t219 * t8 - t39 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t137 * t151, -t97 ^ 2 - t98 ^ 2, t25 * t98 - t26 * t97 + t45, 0, 0, 0, 0, 0, 0, -t50 * t121 + t23, -t22 + t227, -t230 - t231, -t153 * t50 + t47 * t8 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t230 - t231, -t22 - t227, -t221, -t101 * t168 + (-qJD(5) - t121) * t50, t167, -t8 * t121 - t39 * t50 + t2, t121 * t153 + t39 * t47 - t1, 0, 0;];
tauc_reg = t5;
