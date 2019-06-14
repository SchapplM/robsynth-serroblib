% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 01:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:18:32
% EndTime: 2019-05-05 01:18:41
% DurationCPUTime: 3.68s
% Computational Cost: add. (11401->330), mult. (21960->472), div. (0->0), fcn. (15721->12), ass. (0->204)
t186 = sin(qJ(5));
t172 = qJDD(4) + qJDD(5);
t187 = sin(qJ(4));
t190 = cos(qJ(5));
t191 = cos(qJ(4));
t225 = t186 * t191;
t140 = (t187 * t190 + t225) * qJD(2);
t220 = qJD(2) * t191;
t142 = -t186 * t187 * qJD(2) + t190 * t220;
t233 = t142 * t140;
t253 = t172 - t233;
t256 = t186 * t253;
t255 = t190 * t253;
t217 = qJD(2) * qJD(4);
t208 = t191 * t217;
t215 = t187 * qJDD(2);
t149 = -t208 - t215;
t169 = t191 * qJDD(2);
t207 = t187 * t217;
t150 = t169 - t207;
t103 = -t140 * qJD(5) + t186 * t149 + t190 * t150;
t173 = qJD(4) + qJD(5);
t231 = t173 * t140;
t254 = t103 - t231;
t185 = sin(qJ(6));
t204 = -t190 * t149 + t150 * t186;
t102 = -qJD(5) * t142 - t204;
t100 = qJDD(6) - t102;
t189 = cos(qJ(6));
t120 = t142 * t185 - t189 * t173;
t122 = t142 * t189 + t173 * t185;
t98 = t122 * t120;
t249 = t100 - t98;
t252 = t185 * t249;
t251 = t189 * t249;
t179 = sin(pkin(11));
t181 = cos(pkin(11));
t154 = g(1) * t179 - g(2) * t181;
t177 = -g(3) + qJDD(1);
t180 = sin(pkin(6));
t182 = cos(pkin(6));
t250 = t154 * t182 + t177 * t180;
t115 = pkin(5) * t140 - pkin(10) * t142;
t247 = t173 ^ 2;
t160 = qJD(4) * pkin(4) - pkin(9) * t220;
t175 = t187 ^ 2;
t194 = qJD(2) ^ 2;
t228 = t175 * t194;
t127 = -t154 * t180 + t177 * t182;
t178 = qJDD(2) * pkin(2);
t155 = -g(1) * t181 - g(2) * t179;
t188 = sin(qJ(2));
t192 = cos(qJ(2));
t109 = -t155 * t188 + t250 * t192;
t200 = qJDD(3) - t109;
t105 = -t194 * qJ(3) - t178 + t200;
t197 = -qJDD(2) * pkin(8) + t105;
t86 = t191 * t127 + t187 * t197;
t72 = -pkin(4) * t228 + t149 * pkin(9) - qJD(4) * t160 + t86;
t240 = t190 * t72;
t196 = t191 * t197;
t222 = t191 * t194;
t244 = t150 * pkin(9);
t248 = qJDD(4) * pkin(4) + t196 + (-pkin(4) * t222 - pkin(9) * t217 - t127) * t187 - t244;
t42 = t248 * t186 + t240;
t33 = -t247 * pkin(5) + t172 * pkin(10) - t140 * t115 + t42;
t110 = t192 * t155 + t250 * t188;
t203 = 0.2e1 * qJD(3) * qJD(2) + t110;
t216 = qJDD(2) * qJ(3);
t101 = -t194 * pkin(2) + t203 + t216;
t99 = -t194 * pkin(8) + t101;
t80 = -t149 * pkin(4) - pkin(9) * t228 + t160 * t220 + t99;
t43 = -t254 * pkin(10) + (t142 * t173 - t102) * pkin(5) + t80;
t16 = t185 * t33 - t189 * t43;
t17 = t185 * t43 + t189 * t33;
t7 = t185 * t16 + t189 * t17;
t133 = qJD(6) + t140;
t205 = t103 * t185 - t189 * t172;
t65 = (qJD(6) - t133) * t122 + t205;
t118 = t120 ^ 2;
t119 = t122 ^ 2;
t132 = t133 ^ 2;
t138 = t140 ^ 2;
t139 = t142 ^ 2;
t246 = -pkin(8) - pkin(2);
t41 = t186 * t72 - t190 * t248;
t32 = -t172 * pkin(5) - t247 * pkin(10) + t115 * t142 + t41;
t245 = -pkin(5) * t32 + pkin(10) * t7;
t29 = t185 * t32;
t74 = t100 + t98;
t243 = t185 * t74;
t242 = t186 * t80;
t30 = t189 * t32;
t241 = t189 * t74;
t239 = t190 * t80;
t20 = t186 * t42 - t190 * t41;
t238 = t191 * t20;
t113 = t233 + t172;
t237 = t113 * t186;
t236 = t113 * t190;
t235 = t133 * t185;
t234 = t133 * t189;
t230 = t173 * t186;
t229 = t173 * t190;
t176 = t191 ^ 2;
t227 = t176 * t194;
t210 = t187 * t222;
t156 = qJDD(4) + t210;
t224 = t187 * t156;
t157 = qJDD(4) - t210;
t223 = t191 * t157;
t221 = t175 + t176;
t218 = qJD(6) + t133;
t96 = -t119 - t132;
t49 = -t185 * t96 - t241;
t202 = -t103 * t189 - t172 * t185;
t70 = t218 * t120 + t202;
t214 = pkin(5) * t70 + pkin(10) * t49 + t29;
t87 = -t132 - t118;
t46 = t189 * t87 - t252;
t67 = -t218 * t122 - t205;
t213 = pkin(5) * t67 + pkin(10) * t46 - t30;
t212 = t186 * t98;
t211 = t190 * t98;
t209 = -pkin(5) * t190 - pkin(4);
t21 = t186 * t41 + t190 * t42;
t108 = t133 * t120;
t79 = -qJD(6) * t120 - t202;
t69 = t108 + t79;
t37 = t185 * t69 - t189 * t65;
t82 = t118 + t119;
t206 = pkin(5) * t82 + pkin(10) * t37 + t7;
t2 = t186 * t7 - t190 * t32;
t3 = t186 * t32 + t190 * t7;
t1 = t187 * t3 + t191 * t2;
t6 = -t16 * t189 + t17 * t185;
t85 = t187 * t127 - t196;
t51 = t187 * t86 - t191 * t85;
t199 = (-qJD(5) + t173) * t142 - t204;
t193 = qJD(4) ^ 2;
t162 = -t193 - t227;
t161 = -t193 - t228;
t153 = t221 * t194;
t152 = t221 * qJDD(2);
t151 = t169 - 0.2e1 * t207;
t148 = 0.2e1 * t208 + t215;
t135 = (-qJDD(2) * t192 + t188 * t194) * t180;
t134 = (qJDD(2) * t188 + t192 * t194) * t180;
t129 = -t139 + t247;
t128 = t138 - t247;
t126 = -t139 - t247;
t125 = t162 * t191 - t224;
t124 = t161 * t187 + t223;
t123 = t182 * t127;
t116 = t139 - t138;
t111 = -t247 - t138;
t107 = -t119 + t132;
t106 = t118 - t132;
t104 = -t138 - t139;
t97 = t119 - t118;
t95 = -t126 * t186 - t236;
t94 = t126 * t190 - t237;
t93 = t231 + t103;
t88 = (qJD(5) + t173) * t142 + t204;
t84 = t111 * t190 - t256;
t83 = t111 * t186 + t255;
t78 = -qJD(6) * t122 - t205;
t77 = (-t120 * t189 + t122 * t185) * t133;
t76 = (-t120 * t185 - t122 * t189) * t133;
t68 = -t108 + t79;
t62 = -t122 * t235 + t189 * t79;
t61 = t122 * t234 + t185 * t79;
t60 = t120 * t234 - t185 * t78;
t59 = t120 * t235 + t189 * t78;
t58 = t187 * t95 + t191 * t94;
t57 = t186 * t93 + t190 * t199;
t56 = t186 * t199 - t190 * t93;
t55 = t106 * t189 - t243;
t54 = -t107 * t185 + t251;
t53 = t106 * t185 + t241;
t52 = t107 * t189 + t252;
t50 = t187 * t84 + t191 * t83;
t48 = t189 * t96 - t243;
t45 = t185 * t87 + t251;
t38 = -t185 * t68 + t189 * t67;
t36 = t185 * t67 + t189 * t68;
t35 = -t185 * t65 - t189 * t69;
t28 = t187 * t57 + t191 * t56;
t27 = -t186 * t70 + t190 * t49;
t26 = t186 * t49 + t190 * t70;
t25 = -t186 * t67 + t190 * t46;
t24 = t186 * t46 + t190 * t67;
t23 = -t186 * t82 + t190 * t37;
t22 = t186 * t37 + t190 * t82;
t19 = -pkin(10) * t48 + t30;
t18 = -pkin(10) * t45 + t29;
t13 = t187 * t27 + t191 * t26;
t12 = t187 * t25 + t191 * t24;
t11 = -pkin(5) * t48 + t17;
t10 = -pkin(5) * t45 + t16;
t9 = t187 * t23 + t191 * t22;
t8 = t187 * t21 + t238;
t4 = -pkin(10) * t35 - t6;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t177, 0, 0, 0, 0, 0, 0, -t135, -t134, 0, t123 + (t109 * t192 + t110 * t188) * t180, 0, 0, 0, 0, 0, 0, 0, t135, t134, t123 + (t101 * t188 - t105 * t192) * t180, 0, 0, 0, 0, 0, 0, t182 * (-t157 * t187 + t161 * t191) + (-t124 * t192 + t148 * t188) * t180, t182 * (-t156 * t191 - t162 * t187) + (-t125 * t192 + t151 * t188) * t180, (t152 * t192 - t153 * t188) * t180, t182 * (t187 * t85 + t191 * t86) + (t188 * t99 - t192 * t51) * t180, 0, 0, 0, 0, 0, 0, t182 * (-t187 * t83 + t191 * t84) + (t188 * t88 - t192 * t50) * t180, t182 * (-t187 * t94 + t191 * t95) + (t188 * t254 - t192 * t58) * t180, t182 * (-t187 * t56 + t191 * t57) + (t104 * t188 - t192 * t28) * t180, t182 * (-t187 * t20 + t191 * t21) + (t188 * t80 - t192 * t8) * t180, 0, 0, 0, 0, 0, 0, t182 * (-t187 * t24 + t191 * t25) + (-t12 * t192 + t188 * t45) * t180, t182 * (-t187 * t26 + t191 * t27) + (-t13 * t192 + t188 * t48) * t180, t182 * (-t187 * t22 + t191 * t23) + (t188 * t35 - t192 * t9) * t180, t182 * (-t187 * t2 + t191 * t3) + (-t1 * t192 + t188 * t6) * t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t109, -t110, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, -0.2e1 * t178 + t200, t203 + 0.2e1 * t216, -pkin(2) * t105 + qJ(3) * t101, (t150 - t207) * t191, -t148 * t191 - t151 * t187, t223 - t187 * (t193 - t227), (-t149 + t208) * t187, t191 * (-t193 + t228) - t224, 0, qJ(3) * t148 + t246 * t124 + t187 * t99, qJ(3) * t151 + t246 * t125 + t191 * t99, -qJ(3) * t153 - t246 * t152 - t51, qJ(3) * t99 + t246 * t51, t191 * (t103 * t190 - t142 * t230) - t187 * (t103 * t186 + t142 * t229), t191 * (-t186 * t254 - t190 * t88) - t187 * (-t186 * t88 + t190 * t254), t191 * (-t129 * t186 + t255) - t187 * (t129 * t190 + t256), t191 * (-t102 * t186 + t140 * t229) - t187 * (t102 * t190 + t140 * t230), t191 * (t128 * t190 - t237) - t187 * (t128 * t186 + t236), (t191 * (-t140 * t190 + t142 * t186) - t187 * (-t140 * t186 - t142 * t190)) * t173, t191 * (-pkin(9) * t83 + t242) - t187 * (-pkin(4) * t88 + pkin(9) * t84 - t239) + qJ(3) * t88 + t246 * t50, t191 * (-pkin(9) * t94 + t239) - t187 * (-pkin(4) * t254 + pkin(9) * t95 + t242) + qJ(3) * t254 + t246 * t58, t191 * (-pkin(9) * t56 - t20) - t187 * (-pkin(4) * t104 + pkin(9) * t57 + t21) + qJ(3) * t104 + t246 * t28, -pkin(9) * t238 - t187 * (-pkin(4) * t80 + pkin(9) * t21) + qJ(3) * t80 + t246 * t8, t191 * (t190 * t62 + t212) - t187 * (t186 * t62 - t211), t191 * (t186 * t97 + t190 * t38) - t187 * (t186 * t38 - t190 * t97), t191 * (t186 * t69 + t190 * t54) - t187 * (t186 * t54 - t190 * t69), t191 * (t190 * t60 - t212) - t187 * (t186 * t60 + t211), t191 * (-t186 * t65 + t190 * t55) - t187 * (t186 * t55 + t190 * t65), t191 * (t100 * t186 + t190 * t77) - t187 * (-t100 * t190 + t186 * t77), t191 * (-pkin(9) * t24 - t10 * t186 + t18 * t190) - t187 * (-pkin(4) * t45 + pkin(9) * t25 + t10 * t190 + t18 * t186) + qJ(3) * t45 + t246 * t12, t191 * (-pkin(9) * t26 - t11 * t186 + t19 * t190) - t187 * (-pkin(4) * t48 + pkin(9) * t27 + t11 * t190 + t186 * t19) + qJ(3) * t48 + t246 * t13, t191 * (-pkin(9) * t22 + t190 * t4) - t187 * (pkin(9) * t23 + t186 * t4) + t246 * t9 + (pkin(5) * t225 - t187 * t209 + qJ(3)) * t35, (t191 * (pkin(5) * t186 - pkin(10) * t190) - t187 * (-pkin(10) * t186 + t209) + qJ(3)) * t6 + (t246 - pkin(9)) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t194, t105, 0, 0, 0, 0, 0, 0, t124, t125, -t152, t51, 0, 0, 0, 0, 0, 0, t50, t58, t28, t8, 0, 0, 0, 0, 0, 0, t12, t13, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, (-t175 + t176) * t194, t169, -t210, -t215, qJDD(4), -t85, -t86, 0, 0, t233, t116, t93, -t233, t199, t172, pkin(4) * t83 - t41, -t240 - t186 * (-pkin(9) * t207 - t244 - t85) + (-t157 * t186 + t94) * pkin(4), pkin(4) * t56, pkin(4) * t20, t61, t36, t52, t59, t53, t76, pkin(4) * t24 + t213, pkin(4) * t26 + t214, pkin(4) * t22 + t206, pkin(4) * t2 + t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, t116, t93, -t233, t199, t172, -t41, -t42, 0, 0, t61, t36, t52, t59, t53, t76, t213, t214, t206, t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t97, t69, -t98, -t65, t100, -t16, -t17, 0, 0;];
tauJ_reg  = t5;
