% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:17
% EndTime: 2019-03-08 19:16:22
% DurationCPUTime: 2.29s
% Computational Cost: add. (2293->309), mult. (5492->436), div. (0->0), fcn. (4865->16), ass. (0->177)
t148 = cos(qJ(5));
t140 = cos(pkin(12));
t198 = qJD(2) * t140;
t118 = t148 * t198;
t136 = sin(pkin(12));
t145 = sin(qJ(5));
t206 = t145 * t136;
t186 = qJD(2) * t206;
t87 = -t118 + t186;
t86 = qJD(6) + t87;
t236 = t86 - qJD(6);
t137 = sin(pkin(11));
t146 = sin(qJ(2));
t139 = sin(pkin(6));
t200 = qJD(1) * t139;
t188 = t146 * t200;
t109 = t137 * t188;
t141 = cos(pkin(11));
t149 = cos(qJ(2));
t187 = t149 * t200;
t81 = t141 * t187 - t109;
t203 = qJD(4) - t81;
t199 = qJD(2) * t139;
t185 = qJD(1) * t199;
t194 = t139 * qJDD(1);
t235 = t146 * t194 + t149 * t185;
t135 = pkin(12) + qJ(5);
t129 = cos(t135);
t138 = sin(pkin(10));
t142 = cos(pkin(10));
t143 = cos(pkin(6));
t96 = t146 * t137 - t149 * t141;
t158 = t96 * t143;
t167 = t149 * t137 + t146 * t141;
t43 = -t138 * t167 - t142 * t158;
t46 = t138 * t158 - t142 * t167;
t83 = t96 * t139;
t160 = g(1) * t46 + g(2) * t43 - g(3) * t83;
t156 = t160 * t129;
t189 = -t140 * pkin(4) - pkin(3);
t228 = t141 * pkin(2);
t108 = t189 - t228;
t95 = -t148 * t140 + t206;
t97 = t148 * t136 + t145 * t140;
t48 = t95 * pkin(5) - t97 * pkin(9) + t108;
t193 = t140 * qJDD(2);
t195 = t136 * qJDD(2);
t174 = t145 * t195 - t148 * t193;
t91 = t97 * qJD(5);
t59 = qJD(2) * t91 + t174;
t55 = qJDD(6) + t59;
t234 = t48 * t55 - t156;
t115 = t149 * t194;
t82 = qJDD(2) * pkin(2) - t146 * t185 + t115;
t37 = -t137 * t235 + t141 * t82;
t164 = qJDD(4) - t37;
t34 = -qJDD(2) * pkin(3) + t164;
t233 = t160 + t34;
t144 = sin(qJ(6));
t84 = t167 * t139;
t66 = -t84 * t136 + t143 * t140;
t67 = t143 * t136 + t84 * t140;
t22 = t145 * t66 + t148 * t67;
t147 = cos(qJ(6));
t77 = t83 * t147;
t232 = -t144 * t22 + t77;
t89 = t97 * qJD(2);
t121 = t143 * qJD(1) + qJD(3);
t107 = t140 * t121;
t101 = qJD(2) * pkin(2) + t187;
t110 = t141 * t188;
t70 = t137 * t101 + t110;
t68 = qJD(2) * qJ(4) + t70;
t35 = t107 + (-pkin(8) * qJD(2) - t68) * t136;
t40 = t136 * t121 + t140 * t68;
t36 = pkin(8) * t198 + t40;
t168 = t145 * t36 - t148 * t35;
t119 = t143 * qJDD(1) + qJDD(3);
t103 = t140 * t119;
t38 = t137 * t82 + t141 * t235;
t33 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t38;
t14 = t103 + (-pkin(8) * qJDD(2) - t33) * t136;
t18 = t136 * t119 + t140 * t33;
t15 = pkin(8) * t193 + t18;
t170 = t145 * t14 + t148 * t15;
t1 = qJDD(5) * pkin(9) - t168 * qJD(5) + t170;
t69 = t141 * t101 - t109;
t175 = qJD(4) - t69;
t54 = t189 * qJD(2) + t175;
t16 = t87 * pkin(5) - t89 * pkin(9) + t54;
t85 = t167 * t143;
t42 = t138 * t96 - t142 * t85;
t45 = t138 * t85 + t142 * t96;
t161 = g(1) * t45 + g(2) * t42 - g(3) * t84;
t12 = t145 * t35 + t148 * t36;
t171 = -t148 * t14 + t145 * t15;
t2 = -qJDD(5) * pkin(5) + t12 * qJD(5) + t171;
t122 = t137 * pkin(2) + qJ(4);
t222 = pkin(8) + t122;
t92 = t222 * t136;
t93 = t222 * t140;
t50 = t145 * t93 + t148 * t92;
t220 = t50 * qJD(5) + t203 * t95;
t51 = -t145 * t92 + t148 * t93;
t9 = -qJD(5) * pkin(5) + t168;
t90 = t95 * qJD(5);
t231 = -(qJD(6) * t16 + t1) * t95 + t2 * t97 - t9 * t90 + (-qJD(6) * t48 + t220) * t86 - t51 * t55 + t161;
t128 = sin(t135);
t210 = t139 * t142;
t212 = t138 * t139;
t162 = g(1) * (t128 * t45 + t129 * t212) + g(2) * (t128 * t42 - t129 * t210) + g(3) * (-t84 * t128 + t143 * t129);
t230 = (t89 * pkin(5) + pkin(9) * t86) * t86 + t162 + t2;
t190 = qJD(5) * t118 + t145 * t193 + t148 * t195;
t58 = -qJD(5) * t186 + t190;
t73 = t144 * qJD(5) + t147 * t89;
t20 = t73 * qJD(6) - t147 * qJDD(5) + t144 * t58;
t176 = t55 * t97 - t86 * t90;
t197 = qJD(6) * t144;
t192 = t97 * t197;
t229 = -t147 * t176 + t86 * t192;
t196 = t147 * qJD(5);
t71 = t144 * t89 - t196;
t226 = t71 * t86;
t225 = t73 * t86;
t224 = t73 * t89;
t223 = t89 * t71;
t19 = qJD(6) * t196 + t144 * qJDD(5) + t147 * t58 - t89 * t197;
t221 = t19 * t95 + t73 * t91;
t219 = t51 * qJD(5) + t203 * t97;
t78 = t137 * t187 + t110;
t218 = t91 * pkin(5) + t90 * pkin(9) - t78;
t216 = t144 * t55;
t180 = t147 * t86;
t215 = t19 * t144;
t214 = t83 * t144;
t213 = qJD(6) * t97;
t211 = t138 * t146;
t209 = t139 * t149;
t208 = t143 * t146;
t207 = t143 * t149;
t202 = qJDD(1) - g(3);
t201 = t136 ^ 2 + t140 ^ 2;
t191 = t142 * t207;
t177 = -t95 * t20 - t91 * t71;
t10 = qJD(5) * pkin(9) + t12;
t4 = t147 * t10 + t144 * t16;
t173 = t144 * t10 - t147 * t16;
t39 = -t136 * t68 + t107;
t172 = t39 * t136 - t40 * t140;
t169 = t147 * t22 + t214;
t21 = t145 * t67 - t148 * t66;
t79 = qJD(2) * t84;
t166 = qJD(2) * t79 + qJDD(2) * t83;
t165 = t147 * t55 + (-t144 * t87 - t197) * t86;
t163 = -t138 * t207 - t142 * t146;
t159 = -g(1) * t212 + g(2) * t210 - g(3) * t143;
t155 = -pkin(9) * t55 + (-t168 + t9) * t86;
t23 = t189 * qJDD(2) + t164;
t154 = -t176 * t144 - t213 * t180;
t153 = -g(1) * t163 - g(3) * t209;
t125 = -pkin(3) - t228;
t152 = -qJD(2) * t78 + qJDD(2) * t125 + t233;
t150 = qJD(2) ^ 2;
t111 = pkin(2) * t191;
t80 = t96 * t199;
t65 = -qJD(2) * pkin(3) + t175;
t63 = t143 * t128 + t84 * t129;
t61 = -t91 * qJD(5) - t95 * qJDD(5);
t60 = -t90 * qJD(5) + t97 * qJDD(5);
t28 = t128 * t212 - t129 * t45;
t26 = -t128 * t210 - t129 * t42;
t17 = -t136 * t33 + t103;
t8 = t22 * qJD(5) - t97 * t80;
t7 = -t21 * qJD(5) + t95 * t80;
t6 = t59 * pkin(5) - t58 * pkin(9) + t23;
t5 = t147 * t6;
t3 = [t202, 0 (qJDD(2) * t149 - t146 * t150) * t139 (-qJDD(2) * t146 - t149 * t150) * t139, t119 * t143 - t37 * t83 + t38 * t84 - t69 * t79 - t70 * t80 - g(3), -t166 * t140, t166 * t136 (-t136 * t66 + t140 * t67) * qJDD(2) - t201 * t80 * qJD(2), t17 * t66 + t172 * t80 + t18 * t67 + t34 * t83 + t65 * t79 - g(3), 0, 0, 0, 0, 0, -t8 * qJD(5) - t21 * qJDD(5) + t83 * t59 + t79 * t87, -t7 * qJD(5) - t22 * qJDD(5) + t83 * t58 + t79 * t89, 0, 0, 0, 0, 0 (-qJD(6) * t169 - t144 * t7 + t79 * t147) * t86 + t232 * t55 + t8 * t71 + t21 * t20 -(qJD(6) * t232 + t79 * t144 + t147 * t7) * t86 - t169 * t55 + t8 * t73 + t21 * t19; 0, qJDD(2), t115 - g(2) * (t191 - t211) + t153, -g(1) * (t138 * t208 - t142 * t149) - g(2) * (-t138 * t149 - t142 * t208) - t202 * t146 * t139, -g(2) * t111 + t69 * t78 - t70 * t81 + (g(2) * t211 + t38 * t137 + t37 * t141 + t153) * pkin(2), -t152 * t140, t152 * t136, -t17 * t136 + t18 * t140 + t161 + (qJD(2) * t203 + t122 * qJDD(2)) * t201, t34 * t125 - t65 * t78 - g(1) * (pkin(2) * t163 + t46 * pkin(3) - t45 * qJ(4)) - g(2) * (-pkin(2) * t211 + t43 * pkin(3) - t42 * qJ(4) + t111) - g(3) * (pkin(2) * t209 - t83 * pkin(3) + t84 * qJ(4)) + (t18 * t122 + t203 * t40) * t140 + (-t17 * t122 - t203 * t39) * t136, t58 * t97 - t89 * t90, -t58 * t95 - t97 * t59 + t90 * t87 - t89 * t91, t60, t61, 0, -t219 * qJD(5) - t50 * qJDD(5) + t108 * t59 + t23 * t95 + t54 * t91 - t78 * t87 - t156, t220 * qJD(5) - t51 * qJDD(5) + t108 * t58 + t160 * t128 + t23 * t97 - t54 * t90 - t78 * t89, -t73 * t192 + (t19 * t97 - t73 * t90) * t147 -(-t144 * t73 - t147 * t71) * t90 + (-t215 - t147 * t20 + (t144 * t71 - t147 * t73) * qJD(6)) * t97, t221 - t229, t154 + t177, t55 * t95 + t86 * t91, t50 * t20 - t173 * t91 + t5 * t95 + t219 * t71 + (t218 * t86 + (-t10 * t95 - t51 * t86 + t9 * t97) * qJD(6) + t234) * t147 + t231 * t144, t50 * t19 - t4 * t91 + t219 * t73 + (-(-qJD(6) * t10 + t6) * t95 - t9 * t213 + (qJD(6) * t51 - t218) * t86 - t234) * t144 + t231 * t147; 0, 0, 0, 0, t159 + t119, 0, 0, 0, t18 * t136 + t17 * t140 + t159, 0, 0, 0, 0, 0, t61, -t60, 0, 0, 0, 0, 0, t154 - t177, t221 + t229; 0, 0, 0, 0, 0, -t193, t195, -t201 * t150, t172 * qJD(2) + t233, 0, 0, 0, 0, 0, 0.2e1 * t89 * qJD(5) + t174 (-t87 - t186) * qJD(5) + t190, 0, 0, 0, 0, 0, t165 - t223, -t86 ^ 2 * t147 - t216 - t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t87, -t87 ^ 2 + t89 ^ 2 (t87 - t186) * qJD(5) + t190, -t174, qJDD(5), -t54 * t89 - t162 - t171, g(1) * t28 + g(2) * t26 + g(3) * t63 + t54 * t87 - t170, t180 * t73 + t215 (t19 - t226) * t147 + (-t20 - t225) * t144, t180 * t86 + t216 - t224, t165 + t223, -t86 * t89, -pkin(5) * t20 - t12 * t71 + t155 * t144 - t147 * t230 + t173 * t89, -pkin(5) * t19 - t12 * t73 + t144 * t230 + t155 * t147 + t4 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, t19 + t226, -t20 + t225, t55, -t144 * t1 + t5 - t9 * t73 - g(1) * (-t28 * t144 - t46 * t147) - g(2) * (-t26 * t144 - t43 * t147) - g(3) * (-t63 * t144 + t77) + t236 * t4, -t147 * t1 - t144 * t6 + t9 * t71 - g(1) * (t46 * t144 - t28 * t147) - g(2) * (t43 * t144 - t26 * t147) - g(3) * (-t63 * t147 - t214) - t236 * t173;];
tau_reg  = t3;
