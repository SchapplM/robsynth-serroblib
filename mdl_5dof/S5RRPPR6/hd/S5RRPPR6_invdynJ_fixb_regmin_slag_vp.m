% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:47:32
% EndTime: 2021-01-15 19:47:45
% DurationCPUTime: 2.70s
% Computational Cost: add. (3297->364), mult. (7839->495), div. (0->0), fcn. (5805->14), ass. (0->179)
t153 = qJ(2) + pkin(8);
t148 = sin(t153);
t162 = sin(qJ(1));
t165 = cos(qJ(1));
t188 = g(1) * t165 + g(2) * t162;
t150 = cos(t153);
t226 = g(3) * t150;
t173 = t148 * t188 - t226;
t157 = sin(pkin(8));
t215 = cos(pkin(8));
t161 = sin(qJ(2));
t164 = cos(qJ(2));
t159 = -qJ(3) - pkin(6);
t193 = qJD(2) * t159;
t174 = -qJD(3) * t161 + t164 * t193;
t197 = t159 * t161;
t75 = qJDD(2) * pkin(2) + qJD(1) * t174 + qJDD(1) * t197;
t107 = qJD(3) * t164 + t161 * t193;
t134 = t159 * t164;
t84 = qJD(1) * t107 - qJDD(1) * t134;
t33 = -t157 * t84 + t215 * t75;
t32 = -qJDD(2) * pkin(3) + qJDD(4) - t33;
t239 = t32 - t173;
t190 = t215 * t164;
t137 = qJD(1) * t190;
t205 = qJD(1) * t161;
t109 = t157 * t205 - t137;
t106 = qJD(5) + t109;
t160 = sin(qJ(5));
t163 = cos(qJ(5));
t191 = t215 * t161;
t125 = t157 * t164 + t191;
t112 = t125 * qJD(1);
t156 = sin(pkin(9));
t158 = cos(pkin(9));
t96 = qJD(2) * t156 + t112 * t158;
t97 = qJD(2) * t158 - t112 * t156;
t41 = t160 * t96 - t163 * t97;
t238 = t106 * t41;
t126 = t156 * t163 + t158 * t160;
t116 = t126 * qJD(5);
t216 = t126 * t109 + t116;
t180 = -t160 * t97 - t163 * t96;
t237 = t106 * t180;
t235 = g(1) * t162 - g(2) * t165;
t236 = t235 * t148;
t179 = t235 * t150;
t234 = -qJD(5) + t106;
t124 = t156 * t160 - t158 * t163;
t217 = t106 * t124;
t111 = t125 * qJD(2);
t201 = t161 * qJDD(1);
t184 = -qJDD(1) * t190 + t157 * t201;
t81 = qJD(1) * t111 + t184;
t76 = qJDD(5) + t81;
t233 = t106 * t217 - t126 * t76;
t227 = g(3) * t148;
t232 = -t188 * t150 - t227;
t202 = qJD(1) * qJD(2);
t196 = t161 * t202;
t168 = qJDD(1) * t125 - t157 * t196;
t82 = qJD(2) * t137 + t168;
t59 = -qJDD(2) * t158 + t156 * t82;
t60 = qJDD(2) * t156 + t158 * t82;
t9 = -qJD(5) * t180 + t160 * t60 + t163 * t59;
t108 = t109 ^ 2;
t231 = pkin(2) * t161;
t230 = pkin(2) * t164;
t229 = pkin(7) * t158;
t225 = g(3) * t164;
t139 = pkin(2) * t157 + qJ(4);
t224 = pkin(7) + t139;
t144 = pkin(1) + t230;
t105 = pkin(2) * t196 - qJDD(1) * t144 + qJDD(3);
t18 = pkin(3) * t81 - qJ(4) * t82 - qJD(4) * t112 + t105;
t34 = t157 * t75 + t215 * t84;
t29 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t34;
t7 = t156 * t18 + t158 * t29;
t175 = -t157 * t161 + t190;
t114 = t175 * qJD(2);
t223 = qJD(2) * pkin(2);
t199 = t161 * t223;
t44 = pkin(3) * t111 - qJ(4) * t114 - qJD(4) * t125 + t199;
t64 = t107 * t215 + t157 * t174;
t20 = t156 * t44 + t158 * t64;
t132 = -qJD(1) * t144 + qJD(3);
t52 = pkin(3) * t109 - qJ(4) * t112 + t132;
t128 = qJD(1) * t197;
t122 = t128 + t223;
t129 = qJD(1) * t134;
t192 = t215 * t129;
t80 = t122 * t157 - t192;
t72 = qJD(2) * qJ(4) + t80;
t24 = t156 * t52 + t158 * t72;
t62 = pkin(2) * t205 + pkin(3) * t112 + qJ(4) * t109;
t117 = t157 * t129;
t86 = t128 * t215 + t117;
t31 = t156 * t62 + t158 * t86;
t78 = -pkin(3) * t175 - qJ(4) * t125 - t144;
t93 = -t134 * t215 + t157 * t197;
t36 = t156 * t78 + t158 * t93;
t222 = t112 * t41;
t221 = t112 * t180;
t219 = t156 * t81;
t218 = t158 * t81;
t214 = t109 * t156;
t213 = t114 * t156;
t212 = t125 * t156;
t211 = t125 * t158;
t210 = t150 * t162;
t209 = t150 * t165;
t208 = t159 * t162;
t79 = t122 * t215 + t117;
t65 = -qJD(2) * pkin(3) + qJD(4) - t79;
t207 = -qJD(4) + t65;
t154 = t161 ^ 2;
t206 = -t164 ^ 2 + t154;
t204 = qJD(5) * t160;
t203 = qJD(5) * t163;
t200 = t164 * qJDD(1);
t6 = -t156 * t29 + t158 * t18;
t2 = pkin(4) * t81 - pkin(7) * t60 + t6;
t5 = -pkin(7) * t59 + t7;
t198 = -t160 * t5 + t163 * t2;
t19 = -t156 * t64 + t158 * t44;
t23 = -t156 * t72 + t158 * t52;
t30 = -t156 * t86 + t158 * t62;
t35 = -t156 * t93 + t158 * t78;
t63 = t107 * t157 - t174 * t215;
t85 = t128 * t157 - t192;
t92 = -t134 * t157 - t159 * t191;
t189 = -t106 * t216 - t124 * t76;
t142 = -pkin(2) * t215 - pkin(3);
t186 = -t156 * t7 - t158 * t6;
t185 = t160 * t2 + t163 * t5;
t11 = pkin(4) * t109 - pkin(7) * t96 + t23;
t15 = pkin(7) * t97 + t24;
t3 = t11 * t163 - t15 * t160;
t4 = t11 * t160 + t15 * t163;
t183 = -t156 * t23 + t158 * t24;
t21 = -pkin(4) * t175 - pkin(7) * t211 + t35;
t25 = -pkin(7) * t212 + t36;
t182 = -t160 * t25 + t163 * t21;
t181 = t160 * t21 + t163 * t25;
t178 = -0.2e1 * pkin(1) * t202 - pkin(6) * qJDD(2);
t8 = -t160 * t59 + t163 * t60 + t203 * t97 - t204 * t96;
t120 = t224 * t156;
t177 = pkin(7) * t214 - qJD(4) * t158 + qJD(5) * t120 + t31;
t121 = t224 * t158;
t176 = pkin(4) * t112 + qJD(4) * t156 + qJD(5) * t121 + t109 * t229 + t30;
t166 = qJD(2) ^ 2;
t171 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t166 + t235;
t167 = qJD(1) ^ 2;
t170 = pkin(1) * t167 - pkin(6) * qJDD(1) + t188;
t169 = t114 * t65 + t125 * t32 - t188;
t152 = pkin(9) + qJ(5);
t149 = cos(t152);
t147 = sin(t152);
t143 = t159 * t165;
t133 = -pkin(3) * t157 + qJ(4) * t215;
t131 = -pkin(4) * t158 + t142;
t130 = pkin(3) * t215 + qJ(4) * t157 + pkin(2);
t101 = t147 * t162 + t149 * t209;
t100 = -t147 * t209 + t149 * t162;
t99 = t147 * t165 - t149 * t210;
t98 = t147 * t210 + t149 * t165;
t89 = t130 * t164 + t133 * t161 + pkin(1);
t67 = t124 * t125;
t66 = t126 * t125;
t58 = pkin(4) * t212 + t92;
t45 = -pkin(4) * t214 + t85;
t38 = pkin(4) * t213 + t63;
t37 = -pkin(4) * t97 + t65;
t28 = t114 * t126 + t203 * t211 - t204 * t212;
t27 = -t114 * t124 - t116 * t125;
t13 = pkin(4) * t59 + t32;
t12 = -pkin(7) * t213 + t20;
t10 = pkin(4) * t111 - t114 * t229 + t19;
t1 = [qJDD(1), t235, t188, qJDD(1) * t154 + 0.2e1 * t164 * t196, 0.2e1 * t161 * t200 - 0.2e1 * t202 * t206, qJDD(2) * t161 + t164 * t166, qJDD(2) * t164 - t161 * t166, 0, t161 * t178 + t164 * t171, -t161 * t171 + t164 * t178, -qJDD(2) * t92 - t105 * t175 + t111 * t132 - t144 * t81 + t179 + (t109 * t231 - t63) * qJD(2), -qJDD(2) * t93 + t105 * t125 + t114 * t132 - t144 * t82 - t236 + (t112 * t231 - t64) * qJD(2), -t109 * t64 - t111 * t80 + t112 * t63 - t114 * t79 - t125 * t33 + t175 * t34 - t81 * t93 + t82 * t92 - t188, t34 * t93 + t80 * t64 - t33 * t92 - t79 * t63 - t105 * t144 + t132 * t199 - g(1) * (-t144 * t162 - t143) - g(2) * (t144 * t165 - t208), t19 * t109 + t23 * t111 + t156 * t169 + t158 * t179 - t175 * t6 + t35 * t81 + t92 * t59 - t63 * t97, -t20 * t109 - t24 * t111 - t156 * t179 + t158 * t169 + t175 * t7 - t36 * t81 + t92 * t60 + t63 * t96, -t19 * t96 + t20 * t97 - t35 * t60 - t36 * t59 + t236 + t186 * t125 + (-t156 * t24 - t158 * t23) * t114, t7 * t36 + t24 * t20 + t6 * t35 + t23 * t19 + t32 * t92 + t65 * t63 - g(1) * (-t162 * t89 - t143) - g(2) * (t165 * t89 - t208), -t180 * t27 - t67 * t8, t180 * t28 - t27 * t41 - t66 * t8 + t67 * t9, t106 * t27 - t111 * t180 - t175 * t8 - t67 * t76, -t106 * t28 - t111 * t41 + t175 * t9 - t66 * t76, t106 * t111 - t175 * t76, (t10 * t163 - t12 * t160) * t106 + t182 * t76 - t198 * t175 + t3 * t111 + t38 * t41 + t58 * t9 + t13 * t66 + t37 * t28 - g(1) * t99 - g(2) * t101 + (-t106 * t181 + t175 * t4) * qJD(5), -(t10 * t160 + t12 * t163) * t106 - t181 * t76 + t185 * t175 - t4 * t111 - t38 * t180 + t58 * t8 - t13 * t67 + t37 * t27 - g(1) * t98 - g(2) * t100 + (-t106 * t182 + t175 * t3) * qJD(5); 0, 0, 0, -t161 * t167 * t164, t206 * t167, t201, t200, qJDD(2), t161 * t170 - t225, g(3) * t161 + t164 * t170, t85 * qJD(2) - t132 * t112 + (qJDD(2) * t215 - t109 * t205) * pkin(2) + t173 + t33, qJD(2) * t86 + t109 * t132 + (-qJDD(2) * t157 - t112 * t205) * pkin(2) - t34 - t232, (t80 - t85) * t112 + (-t79 + t86) * t109 + (-t157 * t81 - t215 * t82) * pkin(2), t79 * t85 - t80 * t86 + (t215 * t33 - t225 + t157 * t34 + (-qJD(1) * t132 + t188) * t161) * pkin(2), -t139 * t219 - t112 * t23 + t142 * t59 + t85 * t97 + (t156 * t207 - t30) * t109 - t239 * t158, -t139 * t218 + t112 * t24 + t142 * t60 - t85 * t96 + (t158 * t207 + t31) * t109 + t239 * t156, t30 * t96 - t31 * t97 + (qJD(4) * t97 - t109 * t23 - t139 * t59 + t7) * t158 + (qJD(4) * t96 - t109 * t24 + t139 * t60 - t6) * t156 + t232, t32 * t142 - t24 * t31 - t23 * t30 - t65 * t85 - g(3) * (pkin(3) * t150 + qJ(4) * t148 + t230) - t188 * (-t130 * t161 + t133 * t164) + (-t6 * t156 + t7 * t158) * t139 + t183 * qJD(4), t126 * t8 + t180 * t217, -t124 * t8 - t126 * t9 + t180 * t216 + t217 * t41, t221 - t233, t189 + t222, -t106 * t112, (-t120 * t163 - t121 * t160) * t76 + t131 * t9 + t13 * t124 - t3 * t112 - t45 * t41 + t216 * t37 + (t160 * t177 - t163 * t176) * t106 + t173 * t149, -(-t120 * t160 + t121 * t163) * t76 + t131 * t8 + t13 * t126 + t4 * t112 + t45 * t180 - t217 * t37 + (t160 * t176 + t163 * t177) * t106 - t173 * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t112 * qJD(2) + t184, (t137 - t109) * qJD(2) + t168, -t112 ^ 2 - t108, t109 * t80 + t112 * t79 + t105 - t235, -t108 * t156 + t112 * t97 + t218, -t108 * t158 - t112 * t96 - t219, -t156 * t59 - t158 * t60 + (t156 * t96 + t158 * t97) * t109, t109 * t183 - t112 * t65 - t186 - t235, 0, 0, 0, 0, 0, t189 - t222, t221 + t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t96 + t59, t109 * t97 + t60, -t96 ^ 2 - t97 ^ 2, -t125 * t188 + t23 * t96 - t24 * t97 + t226 + t32, 0, 0, 0, 0, 0, t9 - t237, t8 - t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180 * t41, t180 ^ 2 - t41 ^ 2, t8 + t238, -t9 - t237, t76, -g(1) * t100 + g(2) * t98 + t147 * t227 + t180 * t37 + t234 * t4 + t198, g(1) * t101 - g(2) * t99 + t149 * t227 + t234 * t3 + t37 * t41 - t185;];
tau_reg = t1;
