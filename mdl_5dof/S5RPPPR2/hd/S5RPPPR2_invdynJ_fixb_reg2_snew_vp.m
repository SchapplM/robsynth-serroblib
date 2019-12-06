% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:33
% EndTime: 2019-12-05 17:31:46
% DurationCPUTime: 4.58s
% Computational Cost: add. (8467->363), mult. (24370->545), div. (0->0), fcn. (17131->10), ass. (0->212)
t152 = sin(pkin(9));
t153 = sin(pkin(8));
t154 = sin(pkin(7));
t208 = t154 * qJDD(1);
t192 = t153 * t208;
t155 = cos(pkin(9));
t157 = cos(pkin(7));
t156 = cos(pkin(8));
t218 = t154 * t156;
t173 = t152 * t218 + t155 * t157;
t121 = t173 * qJD(1);
t213 = qJD(1) * t154;
t123 = -t152 * t157 * qJD(1) + t155 * t156 * t213;
t224 = t123 * t121;
t248 = t192 - t224;
t251 = t152 * t248;
t250 = t155 * t248;
t149 = t154 ^ 2;
t151 = t157 ^ 2;
t249 = t149 + t151;
t158 = sin(qJ(5));
t242 = t173 * qJDD(1);
t111 = qJDD(5) + t242;
t160 = cos(qJ(5));
t195 = t153 * t213;
t101 = t160 * t123 + t158 * t195;
t99 = t158 * t123 - t160 * t195;
t75 = t101 * t99;
t243 = -t75 + t111;
t247 = t158 * t243;
t246 = t160 * t243;
t162 = qJD(1) ^ 2;
t159 = sin(qJ(1));
t161 = cos(qJ(1));
t178 = t159 * g(2) - t161 * g(3);
t168 = -t162 * pkin(1) + qJDD(1) * qJ(2) + t178;
t240 = 2 * qJD(2);
t245 = qJD(1) * t240 + t168;
t181 = t123 * t195;
t244 = t242 - t181;
t216 = t157 * t162;
t179 = t161 * g(2) + t159 * g(3);
t215 = t162 * qJ(2);
t170 = qJDD(2) - t179 - t215;
t227 = qJDD(1) * pkin(1);
t126 = -t170 + t227;
t241 = t249 * t215 - t126 - t227;
t112 = qJD(5) + t121;
t207 = t157 * qJDD(1);
t140 = t152 * t207;
t209 = qJDD(1) * t156;
t193 = t155 * t209;
t116 = t154 * t193 - t140;
t186 = t158 * t116 - t160 * t192;
t49 = (qJD(5) - t112) * t101 + t186;
t96 = t99 ^ 2;
t97 = t101 ^ 2;
t110 = t112 ^ 2;
t113 = t121 ^ 2;
t114 = t123 ^ 2;
t239 = 2 * qJD(4);
t238 = pkin(4) * t152;
t237 = t153 * pkin(3);
t236 = t157 * g(1);
t177 = -pkin(2) * t157 - qJ(3) * t154;
t174 = -pkin(1) + t177;
t166 = t174 * qJDD(1) + t170;
t180 = -t154 * g(1) + t157 * t245;
t85 = t177 * t216 + t180;
t190 = t153 * t85 - t156 * t166;
t221 = t151 * t162;
t169 = pkin(3) * t207 - qJ(4) * t221 + qJDD(4) + t190;
t176 = -t156 * qJ(4) + t237;
t125 = t176 * t213;
t185 = ((2 * qJD(3)) + t125) * t156;
t45 = t185 * t213 + t169;
t235 = t152 * t45;
t83 = t192 + t224;
t234 = t152 * t83;
t233 = t155 * t45;
t232 = t155 * t83;
t148 = t153 ^ 2;
t222 = t149 * t162;
t142 = t148 * t222;
t167 = (t177 * qJD(1) + t240) * qJD(1);
t187 = qJDD(3) + t236;
t165 = ((-pkin(1) + (-pkin(3) * t156 - qJ(4) * t153) * t157) * t162 + (qJ(2) + t176) * qJDD(1) + t167 + t178) * t154 + t187;
t194 = qJD(3) * t213;
t64 = t156 * t85 + (t166 - 0.2e1 * t194) * t153;
t46 = -pkin(3) * t221 - qJ(4) * t207 - t125 * t195 + t64;
t191 = t152 * t46 - t155 * t165;
t92 = t121 * pkin(4) - t123 * pkin(6);
t21 = -pkin(4) * t192 - pkin(6) * t142 + (t239 + t92) * t123 + t191;
t231 = t158 * t21;
t67 = t75 + t111;
t230 = t158 * t67;
t229 = t160 * t21;
t228 = t160 * t67;
t226 = t112 * t158;
t225 = t112 * t160;
t150 = t156 ^ 2;
t223 = t149 * t150;
t219 = t153 * t162;
t184 = t149 * t156 * t219;
t127 = -t184 + t207;
t220 = t153 * t127;
t128 = -t184 - t207;
t217 = t156 * t128;
t214 = qJD(1) * t153;
t211 = qJD(5) + t112;
t210 = qJDD(1) * t153;
t31 = -t121 * t239 + t152 * t165 + t155 * t46;
t206 = t152 * t75;
t205 = t155 * t75;
t201 = t150 * t222;
t200 = t153 * t224;
t199 = t153 * t216;
t198 = t156 * t216;
t197 = pkin(4) * t155 + pkin(3);
t196 = t121 * t214;
t22 = -pkin(4) * t142 + pkin(6) * t192 - t121 * t92 + t31;
t33 = t242 * pkin(4) - t116 * pkin(6) + (t185 + (pkin(4) * t123 + pkin(6) * t121) * t153) * t213 + t169;
t11 = t158 * t22 - t160 * t33;
t189 = t154 * (t154 * t245 + t236) + t157 * t180;
t183 = t153 * t198;
t182 = t121 * t195;
t12 = t158 * t33 + t160 * t22;
t5 = -t160 * t11 + t158 * t12;
t6 = t158 * t11 + t160 * t12;
t30 = t123 * t239 + t191;
t15 = t152 * t31 - t155 * t30;
t16 = t152 * t30 + t155 * t31;
t63 = 0.2e1 * t156 * t194 + t190;
t34 = t153 * t64 - t156 * t63;
t175 = t149 * t183;
t172 = -t160 * t116 - t158 * t192;
t73 = -t99 * qJD(5) - t172;
t147 = t151 * qJDD(1);
t146 = t149 * qJDD(1);
t141 = t148 * t208;
t133 = t249 * t162;
t131 = (-t151 - t223) * t162;
t130 = -t142 - t221;
t129 = t142 + t201;
t120 = (t199 - t209) * t154;
t119 = (t199 + t209) * t154;
t118 = (-t198 + t210) * t154;
t117 = (t198 + t210) * t154;
t107 = -t114 - t142;
t106 = -t114 + t142;
t105 = t113 - t142;
t95 = t156 * t131 + t220;
t94 = t153 * t130 + t217;
t91 = -t153 * t117 + t156 * t120;
t90 = t140 + (-t193 - t196) * t154;
t89 = -t140 + (t193 - t196) * t154;
t87 = t242 + t181;
t86 = -t142 - t113;
t81 = t187 + (t167 + t168) * t154;
t79 = t112 * t99;
t78 = -t97 + t110;
t77 = t96 - t110;
t76 = -t113 - t114;
t74 = t97 - t96;
t72 = -t101 * qJD(5) - t186;
t71 = -t97 - t110;
t70 = -t152 * t107 - t232;
t69 = t155 * t107 - t234;
t65 = -t110 - t96;
t60 = -t152 * t90 - t155 * t244;
t59 = -t152 * t244 + t155 * t90;
t58 = t96 + t97;
t57 = t155 * t86 - t251;
t56 = t152 * t86 + t250;
t55 = (t101 * t158 - t160 * t99) * t112;
t54 = t211 * t99 + t172;
t53 = t73 + t79;
t52 = t73 - t79;
t50 = -t211 * t101 - t186;
t48 = -t101 * t226 + t160 * t73;
t47 = -t158 * t72 + t99 * t225;
t43 = t153 * t70 - t156 * t89;
t42 = t160 * t77 - t230;
t41 = -t158 * t78 + t246;
t40 = t153 * t57 - t156 * t87;
t39 = t153 * t60 - t156 * t76;
t38 = -t158 * t71 - t228;
t37 = t160 * t71 - t230;
t36 = t160 * t65 - t247;
t35 = t158 * t65 + t246;
t29 = t158 * t53 - t160 * t49;
t28 = -t158 * t52 + t160 * t50;
t27 = -t158 * t49 - t160 * t53;
t26 = -t152 * t54 + t155 * t38;
t25 = t152 * t38 + t155 * t54;
t24 = -t152 * t50 + t155 * t36;
t23 = t152 * t36 + t155 * t50;
t20 = -t152 * t58 + t155 * t29;
t19 = t152 * t29 + t155 * t58;
t18 = t153 * t26 - t156 * t37;
t17 = -pkin(6) * t37 + t229;
t14 = t153 * t24 - t156 * t35;
t13 = -pkin(6) * t35 + t231;
t10 = t153 * t20 - t156 * t27;
t9 = t153 * t16 - t156 * t45;
t8 = -pkin(4) * t37 + t12;
t7 = -pkin(4) * t35 + t11;
t4 = -pkin(6) * t27 - t5;
t3 = t152 * t21 + t155 * t6;
t2 = t152 * t6 - t155 * t21;
t1 = t153 * t3 - t156 * t5;
t32 = [0, 0, 0, 0, 0, qJDD(1), t179, -t178, 0, 0, t146, 0.2e1 * t154 * t207, 0, t147, 0, 0, -t241 * t157, t241 * t154, pkin(1) * t133 + qJ(2) * (t147 + t146) + t189, pkin(1) * t126 + qJ(2) * t189, -t175 + (qJDD(1) * t150 + t183) * t149, t154 * (-t156 * t118 - t153 * t119) + t157 * (t142 - t201), t154 * (t217 - (t151 - t223) * t219) + t157 * t120, t154 * (-t154 * t183 + t141) + t175, t154 * (t156 * (t142 - t221) + t220) + t157 * t117, t147, t154 * (-qJ(3) * t94 + t153 * t81) + t157 * (-pkin(2) * t94 + t63) - pkin(1) * t94 + qJ(2) * (t157 * (-t153 * t128 + t156 * t130) + t154 * t118), t154 * (-qJ(3) * t95 + t156 * t81) + t157 * (-pkin(2) * t95 + t64) - pkin(1) * t95 + qJ(2) * (t157 * (t156 * t127 - t153 * t131) + t154 * t119), -t154 * t34 + qJ(2) * (t157 * (-t156 * t117 - t153 * t120) - t154 * t129) + t174 * t91, qJ(2) * (t157 * (t153 * t63 + t156 * t64) + t154 * t81) + t174 * t34, t154 * (t156 * (t155 * t116 - t152 * t181) + t200) + t157 * (-t152 * t116 - t155 * t181), t154 * (t156 * (-t152 * t89 - t155 * t87) - t153 * (-t114 + t113)) + t157 * (t152 * t87 - t155 * t89), t154 * (t156 * (-t152 * t106 + t250) - t153 * t90) + t157 * (-t155 * t106 - t251), t154 * (t156 * (t152 * t242 + t155 * t182) - t200) + t157 * (-t152 * t182 + t155 * t242), t154 * (t156 * (t155 * t105 - t234) - t153 * t244) + t157 * (-t152 * t105 - t232), (t141 + (t157 * (t121 * t152 + t123 * t155) + (-t121 * t155 + t123 * t152) * t218) * t214) * t154, t154 * (t156 * (-qJ(4) * t56 + t235) - t153 * (-pkin(3) * t56 + t30) - qJ(3) * t40) + t157 * (-pkin(2) * t40 + pkin(3) * t87 - qJ(4) * t57 + t233) - pkin(1) * t40 + qJ(2) * (t157 * (t153 * t87 + t156 * t57) + t154 * t56), t154 * (t156 * (-qJ(4) * t69 + t233) - t153 * (-pkin(3) * t69 + t31) - qJ(3) * t43) + t157 * (-pkin(2) * t43 + pkin(3) * t89 - qJ(4) * t70 - t235) - pkin(1) * t43 + qJ(2) * (t157 * (t153 * t89 + t156 * t70) + t154 * t69), t154 * (t156 * (-qJ(4) * t59 - t15) + t59 * t237 - qJ(3) * t39) + t157 * (-pkin(2) * t39 + pkin(3) * t76 - qJ(4) * t60 - t16) - pkin(1) * t39 + qJ(2) * (t157 * (t153 * t76 + t156 * t60) + t154 * t59), t154 * (-qJ(3) * t9 + t176 * t15) + t157 * (-pkin(2) * t9 + pkin(3) * t45 - qJ(4) * t16) - pkin(1) * t9 + qJ(2) * (t157 * (t153 * t45 + t156 * t16) + t154 * t15), t154 * (t156 * (t155 * t48 + t206) - t153 * (-t101 * t225 - t158 * t73)) + t157 * (-t152 * t48 + t205), t154 * (t156 * (t152 * t74 + t155 * t28) - t153 * (-t158 * t50 - t160 * t52)) + t157 * (-t152 * t28 + t155 * t74), t154 * (t156 * (t152 * t53 + t155 * t41) - t153 * (-t160 * t78 - t247)) + t157 * (-t152 * t41 + t155 * t53), t154 * (t156 * (t155 * t47 - t206) - t153 * (-t160 * t72 - t99 * t226)) + t157 * (-t152 * t47 - t205), t154 * (t156 * (-t152 * t49 + t155 * t42) - t153 * (-t158 * t77 - t228)) + t157 * (-t152 * t42 - t155 * t49), t154 * (t156 * (t152 * t111 + t155 * t55) - t153 * (t101 * t160 + t158 * t99) * t112) + t157 * (t155 * t111 - t152 * t55), t154 * (t156 * (-qJ(4) * t23 + t155 * t13 - t152 * t7) - t153 * (-pkin(3) * t23 - pkin(4) * t50 - pkin(6) * t36 + t229) - qJ(3) * t14) + t157 * (-pkin(2) * t14 + pkin(3) * t35 - qJ(4) * t24 - t152 * t13 - t155 * t7) - pkin(1) * t14 + qJ(2) * (t157 * (t153 * t35 + t156 * t24) + t154 * t23), t154 * (t156 * (-qJ(4) * t25 - t152 * t8 + t155 * t17) - t153 * (-pkin(3) * t25 - pkin(4) * t54 - pkin(6) * t38 - t231) - qJ(3) * t18) + t157 * (-pkin(2) * t18 + pkin(3) * t37 - qJ(4) * t26 - t152 * t17 - t155 * t8) - pkin(1) * t18 + qJ(2) * (t157 * (t153 * t37 + t156 * t26) + t154 * t25), t154 * (t156 * (-qJ(4) * t19 + t155 * t4 + t27 * t238) - t153 * (-pkin(3) * t19 - pkin(4) * t58 - pkin(6) * t29 - t6) - qJ(3) * t10) + t157 * (-pkin(2) * t10 - qJ(4) * t20 - t152 * t4 + t197 * t27) - pkin(1) * t10 + qJ(2) * (t157 * (t153 * t27 + t156 * t20) + t154 * t19), t154 * (t156 * (-qJ(4) * t2 + (-pkin(6) * t155 + t238) * t5) - t153 * (-pkin(3) * t2 + pkin(4) * t21 - pkin(6) * t6) - qJ(3) * t1) + t157 * (-pkin(2) * t1 - qJ(4) * t3 + (pkin(6) * t152 + t197) * t5) - pkin(1) * t1 + qJ(2) * (t157 * (t153 * t5 + t156 * t3) + t154 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, t208, -t133, -t126, 0, 0, 0, 0, 0, 0, t94, t95, t91, t34, 0, 0, 0, 0, 0, 0, t40, t43, t39, t9, 0, 0, 0, 0, 0, 0, t14, t18, t10, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t119, -t129, t81, 0, 0, 0, 0, 0, 0, t56, t69, t59, t15, 0, 0, 0, 0, 0, 0, t23, t25, t19, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t89, t76, t45, 0, 0, 0, 0, 0, 0, t35, t37, t27, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, t53, -t75, -t49, t111, -t11, -t12, 0, 0;];
tauJ_reg = t32;
