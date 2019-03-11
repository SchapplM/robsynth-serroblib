% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:28
% EndTime: 2019-03-09 02:45:34
% DurationCPUTime: 2.32s
% Computational Cost: add. (1615->338), mult. (3019->401), div. (0->0), fcn. (1721->10), ass. (0->190)
t127 = cos(pkin(9));
t84 = -t127 * pkin(1) - pkin(2);
t63 = qJD(1) * t84;
t61 = qJDD(1) * t84;
t131 = sin(qJ(3));
t134 = cos(qJ(3));
t126 = sin(pkin(9));
t83 = t126 * pkin(1) + pkin(7);
t62 = t83 * qJD(1);
t226 = -t134 * qJD(2) + t131 * t62;
t260 = qJD(4) + t226;
t60 = t83 * qJDD(1);
t259 = qJD(2) * qJD(3) + t60;
t248 = pkin(3) + pkin(4);
t258 = t131 * t248;
t257 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t27 = -qJD(3) * pkin(3) + t260;
t203 = qJD(1) * qJD(3);
t187 = t134 * t203;
t200 = t131 * qJDD(1);
t256 = t187 + t200;
t188 = t131 * t203;
t199 = t134 * qJDD(1);
t255 = -t188 + t199;
t193 = t248 * qJD(3);
t205 = t131 * qJD(1);
t25 = -qJ(5) * t205 + t226;
t216 = -qJD(4) - t25;
t19 = -t193 - t216;
t211 = qJD(1) * t134;
t189 = t248 * qJDD(3);
t79 = qJD(6) + t205;
t254 = t79 ^ 2;
t228 = -qJ(5) + t83;
t122 = qJD(3) * qJ(4);
t40 = t131 * qJD(2) + t134 * t62;
t26 = -qJ(5) * t211 + t40;
t23 = -t122 - t26;
t31 = t122 + t40;
t253 = 0.2e1 * t257;
t109 = t131 * qJ(4);
t214 = t134 * pkin(3) + t109;
t118 = qJ(1) + pkin(9);
t100 = sin(t118);
t101 = cos(t118);
t170 = g(1) * t101 + g(2) * t100;
t121 = -pkin(8) - t248;
t171 = t131 * pkin(5) + t134 * pkin(8);
t32 = -pkin(3) * t211 - qJ(4) * t205 + t63;
t22 = pkin(4) * t211 + qJD(5) - t32;
t14 = qJD(1) * t171 + t22;
t209 = qJD(3) * t134;
t179 = -t134 * qJDD(2) + t259 * t131 + t62 * t209;
t167 = -qJDD(4) - t179;
t202 = qJD(1) * qJD(5);
t140 = -qJ(5) * t256 - t131 * t202 - t167;
t186 = -qJD(6) * t14 - t121 * qJDD(3) - t140;
t20 = qJD(3) * pkin(5) - t23;
t114 = t134 * pkin(4);
t44 = -t214 + t84;
t38 = t114 - t44;
t24 = t171 + t38;
t30 = -t131 * qJD(5) + t228 * t209;
t45 = t228 * t131;
t50 = qJDD(6) + t256;
t210 = qJD(3) * t131;
t180 = -t131 * qJDD(2) - t259 * t134 + t62 * t210;
t10 = -t180 + t257;
t204 = qJ(5) * qJDD(1);
t77 = qJ(5) * t188;
t8 = t134 * (t202 + t204) - t10 - t77;
t6 = qJDD(3) * pkin(5) - t8;
t252 = -(qJD(6) * t24 + t30) * t79 + (t20 * qJD(3) + t186) * t131 - t45 * t50 - t6 * t134;
t133 = cos(qJ(6));
t130 = sin(qJ(6));
t207 = qJD(6) * t134;
t192 = t130 * t207;
t233 = t134 * t50;
t251 = -t133 * (t79 * t210 - t233) - t79 * t192;
t206 = qJDD(3) * t83;
t250 = (qJD(1) * t44 + t32) * qJD(3) - t206;
t154 = pkin(5) * t134 + t121 * t131;
t243 = g(3) * t131;
t93 = qJ(4) * t211;
t249 = t134 * t170 + (qJD(1) * t154 + qJD(6) * t121 + t93) * t79 - t6 + t243;
t247 = g(1) * t100;
t244 = g(2) * t101;
t117 = g(3) * t134;
t242 = t24 * t50;
t191 = t130 * t211;
t55 = -t133 * qJD(3) + t191;
t241 = t55 * t79;
t56 = t130 * qJD(3) + t133 * t211;
t240 = t56 * t79;
t15 = qJD(6) * t191 - t130 * qJDD(3) + (-qJD(3) * qJD(6) - t255) * t133;
t239 = t15 * t131 - t56 * t209;
t208 = qJD(6) * t133;
t196 = t79 * t208;
t238 = t130 * t233 + t134 * t196;
t237 = t130 * t50;
t236 = t133 * t50;
t235 = t133 * t56;
t234 = t134 * t23;
t232 = t134 * t55;
t231 = t134 * t56;
t230 = t15 * t130;
t16 = qJD(6) * t56 - t133 * qJDD(3) + t255 * t130;
t229 = t16 * t131;
t105 = t131 * qJD(4);
t227 = qJ(4) * t209 + t105;
t225 = qJ(4) * t134;
t224 = qJDD(3) * pkin(3);
t223 = t100 * t131;
t222 = t100 * t134;
t221 = t101 * t131;
t220 = t101 * t134;
t219 = t130 * t131;
t218 = t131 * t133;
t217 = t23 * qJD(3);
t215 = -qJD(5) - t22;
t123 = t131 ^ 2;
t124 = t134 ^ 2;
t213 = t123 - t124;
t212 = t123 + t124;
t198 = t79 * t219;
t197 = t79 * t218;
t138 = qJD(1) ^ 2;
t195 = t131 * t138 * t134;
t132 = sin(qJ(1));
t190 = -t132 * pkin(1) + t101 * pkin(7);
t17 = t121 * qJD(3) - t216;
t150 = t154 * qJD(3);
t166 = pkin(3) * t199 + t256 * qJ(4) + qJD(1) * t105 - t61;
t152 = pkin(4) * t199 + qJDD(5) + t166;
t2 = qJD(1) * t150 + qJDD(1) * t171 + t152;
t185 = qJD(6) * t17 - t2;
t184 = qJD(1) * t38 + t22;
t181 = t215 * t131;
t176 = 0.2e1 * t187;
t173 = t131 * t193;
t135 = cos(qJ(1));
t172 = t135 * pkin(1) + pkin(3) * t220 + t100 * pkin(7) + (pkin(2) + t109) * t101;
t169 = g(1) * t132 - g(2) * t135;
t137 = qJD(3) ^ 2;
t168 = t137 * t83 + t244;
t165 = -pkin(2) - t214;
t164 = t186 - t117;
t162 = -t198 - t232;
t161 = t198 - t232;
t159 = -t196 - t237;
t158 = -qJD(6) * t130 * t79 + t236;
t157 = -qJD(3) * t226 + t180;
t156 = -g(1) * t221 - g(2) * t223 + t117 + t179;
t151 = -qJDD(4) - t156;
t28 = -t173 + t227;
t9 = -qJD(1) * t173 + t152;
t149 = -qJD(1) * t28 - qJDD(1) * t38 + t244 - t9;
t148 = t168 + 0.2e1 * t61;
t147 = 0.2e1 * t63 * qJD(3) - t206;
t146 = t40 * qJD(3) - t156;
t145 = -t121 * t50 + (-t20 + t26) * t79;
t12 = pkin(3) * t188 - t166;
t43 = pkin(3) * t210 - t227;
t144 = -qJD(1) * t43 - qJDD(1) * t44 - t12 - t168;
t143 = -qJ(5) * t200 - t151;
t11 = -t167 - t224;
t142 = t10 * t134 + t11 * t131 + (-t131 * t31 + t134 * t27) * qJD(3);
t129 = qJ(4) + pkin(5);
t76 = -t123 * t138 - t137;
t71 = g(1) * t222;
t70 = g(1) * t223;
t67 = qJ(4) * t220;
t65 = qJ(4) * t222;
t64 = qJDD(3) + t195;
t59 = qJDD(3) * t134 - t137 * t131;
t58 = qJDD(3) * t131 + t137 * t134;
t57 = pkin(3) * t205 - t93;
t46 = t228 * t134;
t41 = -t248 * t205 + t93;
t36 = -t100 * t130 + t101 * t218;
t35 = -t100 * t133 - t101 * t219;
t34 = -t100 * t218 - t101 * t130;
t33 = t100 * t219 - t101 * t133;
t29 = t134 * qJD(5) + t228 * t210;
t18 = t150 + t227;
t7 = -t189 + t140;
t4 = t130 * t14 + t133 * t17;
t3 = -t130 * t17 + t133 * t14;
t1 = t133 * t2;
t5 = [qJDD(1), t169, g(1) * t135 + g(2) * t132 (t169 + (t126 ^ 2 + t127 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t123 * qJDD(1) + t131 * t176, 0.2e1 * t131 * t199 - 0.2e1 * t213 * t203, t58, t59, 0, t131 * t147 - t134 * t148 + t71, t131 * t148 + t134 * t147 - t70, t131 * t250 + t144 * t134 + t71, t212 * t60 + t142 - t170, t144 * t131 - t134 * t250 + t70, -g(1) * t190 - g(2) * t172 + t12 * t44 + t142 * t83 - t165 * t247 + t32 * t43, t46 * qJDD(3) + t70 + (t134 * t184 - t29) * qJD(3) - t149 * t131, t45 * qJDD(3) - t71 + (t131 * t184 + t30) * qJD(3) + t149 * t134 (-qJD(3) * t19 - qJDD(1) * t46 + t8 + (-qJD(3) * t45 + t29) * qJD(1)) * t134 + (-t217 - qJDD(1) * t45 - t7 + (qJD(3) * t46 - t30) * qJD(1)) * t131 + t170, t7 * t45 + t19 * t30 - t8 * t46 + t23 * t29 + t9 * t38 + t22 * t28 - g(1) * (-t101 * qJ(5) + t190) - g(2) * (pkin(4) * t220 + t172) + (-g(1) * (t165 - t114) + g(2) * qJ(5)) * t100, -t56 * t192 + (-t134 * t15 - t210 * t56) * t133 (t130 * t56 + t133 * t55) * t210 + (t230 - t133 * t16 + (t130 * t55 - t235) * qJD(6)) * t134, t239 - t251, -qJD(3) * t161 + t229 + t238, t50 * t131 + t209 * t79, t3 * t209 - g(1) * t34 - g(2) * t36 + t1 * t131 - t46 * t16 + t29 * t55 + (t18 * t79 + t242 + (-t17 * t131 - t20 * t134 - t45 * t79) * qJD(6)) * t133 + t252 * t130, -t4 * t209 - g(1) * t33 - g(2) * t35 + t46 * t15 + t29 * t56 + (-(-qJD(6) * t45 + t18) * t79 - t242 + t185 * t131 + t20 * t207) * t130 + t252 * t133; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, t59, -t58, t59, 0, t58, t10 * t131 - t11 * t134 - g(3) + (t131 * t27 + t134 * t31) * qJD(3), t58, -t59, 0, -t8 * t131 - t7 * t134 - g(3) + (t131 * t19 - t234) * qJD(3), 0, 0, 0, 0, 0, qJD(3) * t162 - t229 + t238, t239 + t251; 0, 0, 0, 0, -t195, t213 * t138, t200, t199, qJDD(3), -t63 * t205 + t146, t243 + (-qJD(1) * t63 + t170) * t134 + t157, 0.2e1 * t224 - qJDD(4) + (-t131 * t32 + t134 * t57) * qJD(1) + t146 (-pkin(3) * t131 + t225) * qJDD(1) (qJD(1) * t57 - g(3)) * t131 + (qJD(1) * t32 - t170) * t134 - t157 + t253, t10 * qJ(4) - t11 * pkin(3) - t32 * t57 - t27 * t40 - g(1) * (-pkin(3) * t221 + t67) - g(2) * (-pkin(3) * t223 + t65) - g(3) * t214 + t260 * t31, t25 * qJD(3) + t77 + (-qJD(1) * t41 - g(3)) * t131 + (qJD(1) * t215 - t170 - t204) * t134 - t180 + t253, -t26 * qJD(3) - 0.2e1 * t189 + ((-qJ(5) * qJD(3) + t41) * t134 + t181) * qJD(1) + t143 (-t225 + t258) * qJDD(1), -t7 * t248 - t8 * qJ(4) - t19 * t26 - t22 * t41 - g(1) * t67 - g(2) * t65 - g(3) * (t114 + t214) + t216 * t23 + t170 * t258, t235 * t79 - t230 (-t15 - t241) * t133 + (-t16 - t240) * t130 (-t197 + t231) * qJD(1) + t159, qJD(1) * t161 - t158, -t79 * t211, -t129 * t16 + t145 * t130 - t133 * t249 - t3 * t211 + t216 * t55, t129 * t15 + t130 * t249 + t145 * t133 + t4 * t211 + t216 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t200, t76, -t31 * qJD(3) + t205 * t32 - t151 - t224, t76, t64, -t200, t217 - t189 + (-qJ(5) * t209 + t181) * qJD(1) + t143, 0, 0, 0, 0, 0, qJD(3) * t55 - t133 * t254 - t237, qJD(3) * t56 + t130 * t254 - t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176 + t200, 0.2e1 * t188 - t199, -t212 * t138, t247 - t244 + (-t234 + (t19 - t193) * t131) * qJD(1) + t152, 0, 0, 0, 0, 0, qJD(1) * t162 + t158 (-t197 - t231) * qJD(1) + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t55, -t55 ^ 2 + t56 ^ 2, t15 - t241, t16 - t240, t50, -g(1) * t35 + g(2) * t33 + t130 * t164 - t17 * t208 + t20 * t56 + t4 * t79 + t1, g(1) * t36 - g(2) * t34 + t130 * t185 + t133 * t164 - t20 * t55 + t3 * t79;];
tau_reg  = t5;
