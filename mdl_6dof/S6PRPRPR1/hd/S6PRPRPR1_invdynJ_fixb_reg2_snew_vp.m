% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 22:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:11:00
% EndTime: 2019-05-04 22:11:07
% DurationCPUTime: 3.02s
% Computational Cost: add. (13378->358), mult. (27266->539), div. (0->0), fcn. (20135->14), ass. (0->212)
t181 = sin(pkin(12));
t184 = cos(pkin(12));
t193 = cos(qJ(4));
t222 = qJD(2) * t193;
t190 = sin(qJ(4));
t223 = qJD(2) * t190;
t146 = t181 * t223 - t184 * t222;
t148 = t181 * t222 + t184 * t223;
t123 = t148 * t146;
t244 = qJDD(4) - t123;
t249 = t181 * t244;
t248 = t184 * t244;
t189 = sin(qJ(6));
t192 = cos(qJ(6));
t129 = -t192 * qJD(4) + t189 * t148;
t131 = t189 * qJD(4) + t192 * t148;
t101 = t131 * t129;
t216 = qJD(2) * qJD(4);
t210 = t193 * t216;
t215 = t190 * qJDD(2);
t152 = t210 + t215;
t172 = t193 * qJDD(2);
t211 = t190 * t216;
t204 = t172 - t211;
t206 = t181 * t152 - t184 * t204;
t122 = qJDD(6) + t206;
t245 = -t101 + t122;
t247 = t189 * t245;
t246 = t192 * t245;
t219 = t148 * qJD(4);
t102 = t206 + t219;
t232 = sin(pkin(10));
t233 = cos(pkin(10));
t158 = t232 * g(1) - t233 * g(2);
t178 = -g(3) + qJDD(1);
t186 = cos(pkin(6));
t169 = t186 * t178;
t183 = sin(pkin(6));
t136 = -t183 * t158 + qJDD(3) + t169;
t196 = qJD(2) ^ 2;
t191 = sin(qJ(2));
t194 = cos(qJ(2));
t200 = -t233 * g(1) - t232 * g(2);
t227 = t186 * t158;
t203 = t183 * t178 + t227;
t116 = t191 * t203 + t194 * t200;
t114 = -t196 * pkin(2) + t116;
t182 = sin(pkin(11));
t185 = cos(pkin(11));
t115 = -t191 * t200 + t194 * t203;
t198 = qJDD(2) * pkin(2) + t115;
t85 = t185 * t114 + t182 * t198;
t201 = -t196 * pkin(3) + qJDD(2) * pkin(8) + t85;
t199 = t190 * t201;
t225 = t190 * t196;
t197 = -t199 - t152 * qJ(5) + qJDD(4) * pkin(4) + (pkin(4) * t225 + qJ(5) * t216 + t136) * t193;
t159 = qJD(4) * pkin(4) - qJ(5) * t223;
t177 = t193 ^ 2;
t175 = t177 * t196;
t72 = t190 * t136 + t193 * t201;
t62 = -pkin(4) * t175 + t204 * qJ(5) - qJD(4) * t159 + t72;
t35 = -0.2e1 * qJD(5) * t146 + t181 * t197 + t184 * t62;
t143 = qJD(6) + t146;
t125 = t184 * t152 + t181 * t204;
t207 = -t192 * qJDD(4) + t189 * t125;
t75 = (qJD(6) - t143) * t131 + t207;
t127 = t129 ^ 2;
t128 = t131 ^ 2;
t142 = t143 ^ 2;
t144 = t146 ^ 2;
t145 = t148 ^ 2;
t243 = 0.2e1 * qJD(5);
t242 = pkin(5) * t181;
t208 = t182 * t114 - t185 * t198;
t83 = -qJDD(2) * pkin(3) - t196 * pkin(8) + t208;
t67 = -t204 * pkin(4) - qJ(5) * t175 + t159 * t223 + qJDD(5) + t83;
t240 = t181 * t67;
t239 = t184 * t67;
t117 = t146 * pkin(5) - t148 * pkin(9);
t195 = qJD(4) ^ 2;
t209 = t181 * t62 - t184 * t197;
t30 = -qJDD(4) * pkin(5) - t195 * pkin(9) + (t243 + t117) * t148 + t209;
t238 = t189 * t30;
t88 = t101 + t122;
t237 = t189 * t88;
t34 = t148 * t243 + t209;
t17 = t181 * t35 - t184 * t34;
t236 = t190 * t17;
t235 = t192 * t30;
t234 = t192 * t88;
t231 = t143 * t189;
t230 = t143 * t192;
t120 = qJDD(4) + t123;
t229 = t181 * t120;
t228 = t184 * t120;
t167 = t193 * t225;
t160 = qJDD(4) + t167;
t226 = t190 * t160;
t161 = qJDD(4) - t167;
t224 = t193 * t161;
t220 = t146 * qJD(4);
t217 = qJD(6) + t143;
t214 = t181 * t101;
t213 = t184 * t101;
t212 = -pkin(5) * t184 - pkin(4);
t18 = t181 * t34 + t184 * t35;
t31 = -t195 * pkin(5) + qJDD(4) * pkin(9) - t146 * t117 + t35;
t205 = -t125 + t220;
t46 = t102 * pkin(5) + t205 * pkin(9) + t67;
t20 = t189 * t31 - t192 * t46;
t21 = t189 * t46 + t192 * t31;
t11 = t189 * t20 + t192 * t21;
t71 = -t193 * t136 + t199;
t44 = t190 * t71 + t193 * t72;
t4 = t181 * t11 - t184 * t30;
t5 = t184 * t11 + t181 * t30;
t3 = -t190 * t4 + t193 * t5;
t153 = t172 - 0.2e1 * t211;
t10 = t189 * t21 - t192 * t20;
t202 = -t189 * qJDD(4) - t192 * t125;
t103 = -t206 + t219;
t95 = -t129 * qJD(6) - t202;
t176 = t190 ^ 2;
t173 = t176 * t196;
t166 = -t175 - t195;
t165 = -t173 - t195;
t157 = t173 + t175;
t156 = (t176 + t177) * qJDD(2);
t155 = -t182 * qJDD(2) - t185 * t196;
t154 = t185 * qJDD(2) - t182 * t196;
t151 = 0.2e1 * t210 + t215;
t139 = -t145 - t195;
t138 = -t145 + t195;
t137 = t144 - t195;
t135 = -t190 * t165 - t224;
t134 = t193 * t166 - t226;
t133 = -t190 * t161 + t193 * t165;
t132 = t193 * t160 + t190 * t166;
t126 = t182 * t156 + t185 * t157;
t118 = -t195 - t144;
t113 = t143 * t129;
t112 = t182 * t135 - t185 * t151;
t111 = t182 * t134 + t185 * t153;
t110 = -t128 + t142;
t109 = t127 - t142;
t105 = t125 + t220;
t100 = -t144 - t145;
t99 = t128 - t127;
t98 = -t128 - t142;
t97 = -t181 * t139 - t228;
t96 = t184 * t139 - t229;
t94 = -t131 * qJD(6) - t207;
t93 = -t142 - t127;
t92 = t127 + t128;
t91 = t184 * t118 - t249;
t90 = t181 * t118 + t248;
t86 = (-t129 * t192 + t131 * t189) * t143;
t82 = t184 * t103 + t181 * t105;
t81 = t181 * t103 - t184 * t105;
t80 = t217 * t129 + t202;
t79 = t113 + t95;
t78 = -t113 + t95;
t77 = -t217 * t131 - t207;
t74 = -t131 * t231 + t192 * t95;
t73 = t129 * t230 - t189 * t94;
t69 = -t190 * t96 + t193 * t97;
t68 = t190 * t97 + t193 * t96;
t66 = t192 * t109 - t237;
t65 = -t189 * t110 + t246;
t64 = -t189 * t98 - t234;
t63 = t192 * t98 - t237;
t61 = t192 * t93 - t247;
t60 = t189 * t93 + t246;
t56 = -t190 * t90 + t193 * t91;
t55 = t190 * t91 + t193 * t90;
t54 = t182 * t69 + t185 * t205;
t53 = t182 * t85 - t185 * t208;
t52 = -t190 * t81 + t193 * t82;
t51 = t190 * t82 + t193 * t81;
t50 = -t185 * t102 + t182 * t56;
t49 = -t189 * t78 + t192 * t77;
t48 = t189 * t79 - t192 * t75;
t47 = -t189 * t75 - t192 * t79;
t43 = t190 * t72 - t193 * t71;
t42 = -t181 * t80 + t184 * t64;
t41 = t181 * t64 + t184 * t80;
t40 = -t181 * t77 + t184 * t61;
t39 = t181 * t61 + t184 * t77;
t38 = -t185 * t100 + t182 * t52;
t37 = -t181 * t92 + t184 * t48;
t36 = t181 * t48 + t184 * t92;
t32 = t182 * t44 - t185 * t83;
t29 = -t190 * t41 + t193 * t42;
t28 = t190 * t42 + t193 * t41;
t27 = -t190 * t39 + t193 * t40;
t26 = t190 * t40 + t193 * t39;
t25 = -pkin(9) * t63 + t235;
t24 = -pkin(9) * t60 + t238;
t23 = -t190 * t36 + t193 * t37;
t22 = t190 * t37 + t193 * t36;
t16 = t182 * t29 - t185 * t63;
t15 = t182 * t27 - t185 * t60;
t14 = -pkin(5) * t63 + t21;
t13 = -pkin(5) * t60 + t20;
t12 = t182 * t23 - t185 * t47;
t9 = t193 * t18 - t236;
t8 = t193 * t17 + t190 * t18;
t7 = -pkin(9) * t47 - t10;
t6 = t182 * t9 - t185 * t67;
t2 = t190 * t5 + t193 * t4;
t1 = -t185 * t10 + t182 * t3;
t19 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t178, 0, 0, 0, 0, 0, 0, (qJDD(2) * t194 - t191 * t196) * t183, (-qJDD(2) * t191 - t194 * t196) * t183, 0, t186 * t169 + (t194 * t115 + t191 * t116 - t227) * t183, 0, 0, 0, 0, 0, 0, (t154 * t194 + t155 * t191) * t183, (-t154 * t191 + t155 * t194) * t183, 0, t186 * t136 + (t191 * (t182 * t208 + t185 * t85) + t194 * t53) * t183, 0, 0, 0, 0, 0, 0, t186 * t132 + (t191 * (t185 * t134 - t182 * t153) + t194 * t111) * t183, t186 * t133 + (t191 * (t185 * t135 + t182 * t151) + t194 * t112) * t183, (t191 * (t185 * t156 - t182 * t157) + t194 * t126) * t183, t186 * t43 + (t191 * (t182 * t83 + t185 * t44) + t194 * t32) * t183, 0, 0, 0, 0, 0, 0, t186 * t55 + (t191 * (t182 * t102 + t185 * t56) + t194 * t50) * t183, t186 * t68 + (t191 * (-t182 * t205 + t185 * t69) + t194 * t54) * t183, t186 * t51 + (t191 * (t182 * t100 + t185 * t52) + t194 * t38) * t183, t186 * t8 + (t191 * (t182 * t67 + t185 * t9) + t194 * t6) * t183, 0, 0, 0, 0, 0, 0, t186 * t26 + (t191 * (t182 * t60 + t185 * t27) + t194 * t15) * t183, t186 * t28 + (t191 * (t182 * t63 + t185 * t29) + t194 * t16) * t183, t186 * t22 + (t191 * (t182 * t47 + t185 * t23) + t194 * t12) * t183, t186 * t2 + (t191 * (t182 * t10 + t185 * t3) + t194 * t1) * t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t115, -t116, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t154 - t208, pkin(2) * t155 - t85, 0, pkin(2) * t53, (t152 + t210) * t190, t193 * t151 + t190 * t153, t226 + t193 * (-t173 + t195), t153 * t193, t190 * (t175 - t195) + t224, 0, pkin(2) * t111 + pkin(3) * t153 + pkin(8) * t134 - t193 * t83, pkin(2) * t112 - pkin(3) * t151 + pkin(8) * t135 + t190 * t83, pkin(2) * t126 + pkin(3) * t157 + pkin(8) * t156 + t44, pkin(2) * t32 - pkin(3) * t83 + pkin(8) * t44, t190 * (t184 * t125 - t181 * t219) + t193 * (t181 * t125 + t184 * t219), t190 * (-t184 * t102 + t181 * t205) + t193 * (-t181 * t102 - t184 * t205), t190 * (-t181 * t138 + t248) + t193 * (t184 * t138 + t249), t190 * (t181 * t206 + t184 * t220) + t193 * (t181 * t220 - t184 * t206), t190 * (t184 * t137 - t229) + t193 * (t181 * t137 + t228), (t190 * (-t146 * t184 + t148 * t181) + t193 * (-t146 * t181 - t148 * t184)) * qJD(4), t190 * (-qJ(5) * t90 + t240) + t193 * (-pkin(4) * t102 + qJ(5) * t91 - t239) - pkin(3) * t102 + pkin(8) * t56 + pkin(2) * t50, t190 * (-qJ(5) * t96 + t239) + t193 * (pkin(4) * t205 + qJ(5) * t97 + t240) + pkin(3) * t205 + pkin(8) * t69 + pkin(2) * t54, t190 * (-qJ(5) * t81 - t17) + t193 * (-pkin(4) * t100 + qJ(5) * t82 + t18) - pkin(3) * t100 + pkin(8) * t52 + pkin(2) * t38, -qJ(5) * t236 + t193 * (-pkin(4) * t67 + qJ(5) * t18) - pkin(3) * t67 + pkin(8) * t9 + pkin(2) * t6, t190 * (t184 * t74 + t214) + t193 * (t181 * t74 - t213), t190 * (t181 * t99 + t184 * t49) + t193 * (t181 * t49 - t184 * t99), t190 * (t181 * t79 + t184 * t65) + t193 * (t181 * t65 - t184 * t79), t190 * (t184 * t73 - t214) + t193 * (t181 * t73 + t213), t190 * (-t181 * t75 + t184 * t66) + t193 * (t181 * t66 + t184 * t75), t190 * (t181 * t122 + t184 * t86) + t193 * (-t184 * t122 + t181 * t86), t190 * (-qJ(5) * t39 - t181 * t13 + t184 * t24) + t193 * (-pkin(4) * t60 + qJ(5) * t40 + t184 * t13 + t181 * t24) - pkin(3) * t60 + pkin(8) * t27 + pkin(2) * t15, t190 * (-qJ(5) * t41 - t181 * t14 + t184 * t25) + t193 * (-pkin(4) * t63 + qJ(5) * t42 + t184 * t14 + t181 * t25) - pkin(3) * t63 + pkin(8) * t29 + pkin(2) * t16, t190 * (-qJ(5) * t36 + t184 * t7) + t193 * (qJ(5) * t37 + t181 * t7) + pkin(8) * t23 + pkin(2) * t12 + (t190 * t242 + t193 * t212 - pkin(3)) * t47, pkin(2) * t1 + (t190 * (-pkin(9) * t184 + t242) + t193 * (-pkin(9) * t181 + t212) - pkin(3)) * t10 + (pkin(8) + qJ(5)) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, 0, 0, 0, 0, 0, t132, t133, 0, t43, 0, 0, 0, 0, 0, 0, t55, t68, t51, t8, 0, 0, 0, 0, 0, 0, t26, t28, t22, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167, t173 - t175, t215, t167, t172, qJDD(4), -t71, -t72, 0, 0, t123, t145 - t144, t105, -t123, t103, qJDD(4), pkin(4) * t90 - t34, pkin(4) * t96 - t35, pkin(4) * t81, pkin(4) * t17, t131 * t230 + t189 * t95, t189 * t77 + t192 * t78, t192 * t110 + t247, t129 * t231 + t192 * t94, t189 * t109 + t234, (-t129 * t189 - t131 * t192) * t143, pkin(4) * t39 + pkin(5) * t77 + pkin(9) * t61 - t235, pkin(4) * t41 + pkin(5) * t80 + pkin(9) * t64 + t238, pkin(4) * t36 + pkin(5) * t92 + pkin(9) * t48 + t11, pkin(4) * t4 - pkin(5) * t30 + pkin(9) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, -t205, t100, t67, 0, 0, 0, 0, 0, 0, t60, t63, t47, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t99, t79, -t101, -t75, t122, -t20, -t21, 0, 0;];
tauJ_reg  = t19;
