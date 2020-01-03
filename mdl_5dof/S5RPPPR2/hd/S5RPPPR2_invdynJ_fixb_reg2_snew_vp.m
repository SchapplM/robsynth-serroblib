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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:22:57
% EndTime: 2020-01-03 11:23:13
% DurationCPUTime: 4.70s
% Computational Cost: add. (8467->360), mult. (24370->545), div. (0->0), fcn. (17131->10), ass. (0->211)
t154 = sin(pkin(9));
t155 = sin(pkin(8));
t156 = sin(pkin(7));
t209 = t156 * qJDD(1);
t193 = t155 * t209;
t157 = cos(pkin(9));
t159 = cos(pkin(7));
t158 = cos(pkin(8));
t218 = t156 * t158;
t174 = t154 * t218 + t157 * t159;
t121 = t174 * qJD(1);
t214 = qJD(1) * t156;
t123 = -t154 * t159 * qJD(1) + t157 * t158 * t214;
t224 = t123 * t121;
t247 = t193 - t224;
t250 = t154 * t247;
t249 = t157 * t247;
t149 = t156 ^ 2;
t151 = t159 ^ 2;
t248 = t149 + t151;
t160 = sin(qJ(5));
t241 = t174 * qJDD(1);
t111 = qJDD(5) + t241;
t162 = cos(qJ(5));
t196 = t155 * t214;
t101 = t162 * t123 + t160 * t196;
t99 = t160 * t123 - t162 * t196;
t75 = t101 * t99;
t242 = -t75 + t111;
t246 = t160 * t242;
t245 = t162 * t242;
t164 = qJD(1) ^ 2;
t161 = sin(qJ(1));
t163 = cos(qJ(1));
t179 = t161 * g(2) - t163 * g(3);
t170 = -t164 * pkin(1) + qJDD(1) * qJ(2) - t179;
t239 = 2 * qJD(2);
t244 = qJD(1) * t239 + t170;
t182 = t123 * t196;
t243 = t241 - t182;
t216 = t159 * t164;
t152 = t164 * qJ(2);
t153 = qJDD(1) * pkin(1);
t180 = -t163 * g(2) - t161 * g(3);
t126 = qJDD(2) - t152 - t153 - t180;
t240 = t248 * t152 + t126 - t153;
t112 = qJD(5) + t121;
t208 = t159 * qJDD(1);
t140 = t154 * t208;
t210 = qJDD(1) * t158;
t194 = t157 * t210;
t116 = t156 * t194 - t140;
t187 = t160 * t116 - t162 * t193;
t49 = (qJD(5) - t112) * t101 + t187;
t96 = t99 ^ 2;
t97 = t101 ^ 2;
t110 = t112 ^ 2;
t113 = t121 ^ 2;
t114 = t123 ^ 2;
t238 = 2 * qJD(4);
t237 = pkin(4) * t154;
t236 = t155 * pkin(3);
t235 = t159 * g(1);
t178 = -pkin(2) * t159 - qJ(3) * t156;
t168 = qJDD(1) * t178 + t126;
t181 = -t156 * g(1) + t159 * t244;
t85 = t178 * t216 + t181;
t191 = t155 * t85 - t158 * t168;
t221 = t151 * t164;
t171 = pkin(3) * t208 - qJ(4) * t221 + qJDD(4) + t191;
t177 = -t158 * qJ(4) + t236;
t125 = t177 * t214;
t186 = ((2 * qJD(3)) + t125) * t158;
t45 = t186 * t214 + t171;
t234 = t154 * t45;
t83 = t193 + t224;
t233 = t154 * t83;
t232 = t157 * t45;
t231 = t157 * t83;
t148 = t155 ^ 2;
t222 = t149 * t164;
t142 = t148 * t222;
t169 = (t178 * qJD(1) + t239) * qJD(1);
t188 = qJDD(3) + t235;
t167 = ((-pkin(1) + (-pkin(3) * t158 - qJ(4) * t155) * t159) * t164 + (qJ(2) + t177) * qJDD(1) + t169 - t179) * t156 + t188;
t195 = qJD(3) * t214;
t64 = t158 * t85 + (t168 - 0.2e1 * t195) * t155;
t46 = -pkin(3) * t221 - qJ(4) * t208 - t125 * t196 + t64;
t192 = t154 * t46 - t157 * t167;
t92 = t121 * pkin(4) - t123 * pkin(6);
t21 = -pkin(4) * t193 - pkin(6) * t142 + (t238 + t92) * t123 + t192;
t230 = t160 * t21;
t67 = t75 + t111;
t229 = t160 * t67;
t228 = t162 * t21;
t227 = t162 * t67;
t226 = t112 * t160;
t225 = t112 * t162;
t150 = t158 ^ 2;
t223 = t149 * t150;
t219 = t155 * t164;
t185 = t149 * t158 * t219;
t127 = -t185 + t208;
t220 = t155 * t127;
t128 = -t185 - t208;
t217 = t158 * t128;
t215 = qJD(1) * t155;
t212 = qJD(5) + t112;
t211 = qJDD(1) * t155;
t31 = -t121 * t238 + t154 * t167 + t157 * t46;
t207 = t154 * t75;
t206 = t157 * t75;
t202 = t150 * t222;
t201 = t155 * t224;
t200 = t155 * t216;
t199 = t158 * t216;
t198 = pkin(4) * t157 + pkin(3);
t197 = t121 * t215;
t22 = -pkin(4) * t142 + pkin(6) * t193 - t121 * t92 + t31;
t33 = t241 * pkin(4) - t116 * pkin(6) + (t186 + (pkin(4) * t123 + pkin(6) * t121) * t155) * t214 + t171;
t11 = t160 * t22 - t162 * t33;
t190 = t156 * (t156 * t244 + t235) + t159 * t181;
t184 = t155 * t199;
t183 = t121 * t196;
t12 = t160 * t33 + t162 * t22;
t5 = -t162 * t11 + t160 * t12;
t6 = t160 * t11 + t162 * t12;
t30 = t123 * t238 + t192;
t15 = t154 * t31 - t157 * t30;
t16 = t154 * t30 + t157 * t31;
t63 = 0.2e1 * t158 * t195 + t191;
t34 = t155 * t64 - t158 * t63;
t176 = t149 * t184;
t175 = -pkin(1) + t178;
t173 = -t162 * t116 - t160 * t193;
t73 = -t99 * qJD(5) - t173;
t147 = t151 * qJDD(1);
t146 = t149 * qJDD(1);
t141 = t148 * t209;
t133 = t248 * t164;
t131 = (-t151 - t223) * t164;
t130 = -t142 - t221;
t129 = t142 + t202;
t120 = (t200 - t210) * t156;
t119 = (t200 + t210) * t156;
t118 = (-t199 + t211) * t156;
t117 = (t199 + t211) * t156;
t107 = -t114 - t142;
t106 = -t114 + t142;
t105 = t113 - t142;
t95 = t158 * t131 + t220;
t94 = t155 * t130 + t217;
t91 = -t155 * t117 + t158 * t120;
t90 = t140 + (-t194 - t197) * t156;
t89 = -t140 + (t194 - t197) * t156;
t87 = t241 + t182;
t86 = -t142 - t113;
t81 = t188 + (t169 + t170) * t156;
t79 = t112 * t99;
t78 = -t97 + t110;
t77 = t96 - t110;
t76 = -t113 - t114;
t74 = t97 - t96;
t72 = -t101 * qJD(5) - t187;
t71 = -t97 - t110;
t70 = -t154 * t107 - t231;
t69 = t157 * t107 - t233;
t65 = -t110 - t96;
t60 = -t154 * t90 - t157 * t243;
t59 = -t154 * t243 + t157 * t90;
t58 = t96 + t97;
t57 = t157 * t86 - t250;
t56 = t154 * t86 + t249;
t55 = (t101 * t160 - t162 * t99) * t112;
t54 = t212 * t99 + t173;
t53 = t73 + t79;
t52 = t73 - t79;
t50 = -t212 * t101 - t187;
t48 = -t101 * t226 + t162 * t73;
t47 = -t160 * t72 + t99 * t225;
t43 = t155 * t70 - t158 * t89;
t42 = t162 * t77 - t229;
t41 = -t160 * t78 + t245;
t40 = t155 * t57 - t158 * t87;
t39 = t155 * t60 - t158 * t76;
t38 = -t160 * t71 - t227;
t37 = t162 * t71 - t229;
t36 = t162 * t65 - t246;
t35 = t160 * t65 + t245;
t29 = t160 * t53 - t162 * t49;
t28 = -t160 * t52 + t162 * t50;
t27 = -t160 * t49 - t162 * t53;
t26 = -t154 * t54 + t157 * t38;
t25 = t154 * t38 + t157 * t54;
t24 = -t154 * t50 + t157 * t36;
t23 = t154 * t36 + t157 * t50;
t20 = -t154 * t58 + t157 * t29;
t19 = t154 * t29 + t157 * t58;
t18 = t155 * t26 - t158 * t37;
t17 = -pkin(6) * t37 + t228;
t14 = t155 * t24 - t158 * t35;
t13 = -pkin(6) * t35 + t230;
t10 = t155 * t20 - t158 * t27;
t9 = t155 * t16 - t158 * t45;
t8 = -pkin(4) * t37 + t12;
t7 = -pkin(4) * t35 + t11;
t4 = -pkin(6) * t27 - t5;
t3 = t154 * t21 + t157 * t6;
t2 = t154 * t6 - t157 * t21;
t1 = t155 * t3 - t158 * t5;
t32 = [0, 0, 0, 0, 0, qJDD(1), t180, t179, 0, 0, t146, 0.2e1 * t156 * t208, 0, t147, 0, 0, -t240 * t159, t240 * t156, pkin(1) * t133 + qJ(2) * (t147 + t146) + t190, -pkin(1) * t126 + qJ(2) * t190, -t176 + (qJDD(1) * t150 + t184) * t149, t156 * (-t158 * t118 - t155 * t119) + t159 * (t142 - t202), t156 * (t217 - (t151 - t223) * t219) + t159 * t120, t156 * (-t156 * t184 + t141) + t176, t156 * (t158 * (t142 - t221) + t220) + t159 * t117, t147, t156 * (-qJ(3) * t94 + t155 * t81) + t159 * (-pkin(2) * t94 + t63) - pkin(1) * t94 + qJ(2) * (t159 * (-t155 * t128 + t158 * t130) + t156 * t118), t156 * (-qJ(3) * t95 + t158 * t81) + t159 * (-pkin(2) * t95 + t64) - pkin(1) * t95 + qJ(2) * (t159 * (t158 * t127 - t155 * t131) + t156 * t119), -t156 * t34 + qJ(2) * (t159 * (-t158 * t117 - t155 * t120) - t156 * t129) + t175 * t91, qJ(2) * (t159 * (t155 * t63 + t158 * t64) + t156 * t81) + t175 * t34, t156 * (t158 * (t157 * t116 - t154 * t182) + t201) + t159 * (-t154 * t116 - t157 * t182), t156 * (t158 * (-t154 * t89 - t157 * t87) - t155 * (-t114 + t113)) + t159 * (t154 * t87 - t157 * t89), t156 * (t158 * (-t154 * t106 + t249) - t155 * t90) + t159 * (-t157 * t106 - t250), t156 * (t158 * (t154 * t241 + t157 * t183) - t201) + t159 * (-t154 * t183 + t157 * t241), t156 * (t158 * (t157 * t105 - t233) - t155 * t243) + t159 * (-t154 * t105 - t231), (t141 + (t159 * (t121 * t154 + t123 * t157) + (-t121 * t157 + t123 * t154) * t218) * t215) * t156, t156 * (t158 * (-qJ(4) * t56 + t234) - t155 * (-pkin(3) * t56 + t30) - qJ(3) * t40) + t159 * (-pkin(2) * t40 + pkin(3) * t87 - qJ(4) * t57 + t232) - pkin(1) * t40 + qJ(2) * (t159 * (t155 * t87 + t158 * t57) + t156 * t56), t156 * (t158 * (-qJ(4) * t69 + t232) - t155 * (-pkin(3) * t69 + t31) - qJ(3) * t43) + t159 * (-pkin(2) * t43 + pkin(3) * t89 - qJ(4) * t70 - t234) - pkin(1) * t43 + qJ(2) * (t159 * (t155 * t89 + t158 * t70) + t156 * t69), t156 * (t158 * (-qJ(4) * t59 - t15) + t59 * t236 - qJ(3) * t39) + t159 * (-pkin(2) * t39 + pkin(3) * t76 - qJ(4) * t60 - t16) - pkin(1) * t39 + qJ(2) * (t159 * (t155 * t76 + t158 * t60) + t156 * t59), t156 * (-qJ(3) * t9 + t15 * t177) + t159 * (-pkin(2) * t9 + pkin(3) * t45 - qJ(4) * t16) - pkin(1) * t9 + qJ(2) * (t159 * (t155 * t45 + t158 * t16) + t156 * t15), t156 * (t158 * (t157 * t48 + t207) - t155 * (-t101 * t225 - t160 * t73)) + t159 * (-t154 * t48 + t206), t156 * (t158 * (t154 * t74 + t157 * t28) - t155 * (-t160 * t50 - t162 * t52)) + t159 * (-t154 * t28 + t157 * t74), t156 * (t158 * (t154 * t53 + t157 * t41) - t155 * (-t162 * t78 - t246)) + t159 * (-t154 * t41 + t157 * t53), t156 * (t158 * (t157 * t47 - t207) - t155 * (-t162 * t72 - t226 * t99)) + t159 * (-t154 * t47 - t206), t156 * (t158 * (-t154 * t49 + t157 * t42) - t155 * (-t160 * t77 - t227)) + t159 * (-t154 * t42 - t157 * t49), t156 * (t158 * (t154 * t111 + t157 * t55) - t155 * (t101 * t162 + t160 * t99) * t112) + t159 * (t157 * t111 - t154 * t55), t156 * (t158 * (-qJ(4) * t23 + t157 * t13 - t154 * t7) - t155 * (-pkin(3) * t23 - pkin(4) * t50 - pkin(6) * t36 + t228) - qJ(3) * t14) + t159 * (-pkin(2) * t14 + pkin(3) * t35 - qJ(4) * t24 - t154 * t13 - t157 * t7) - pkin(1) * t14 + qJ(2) * (t159 * (t155 * t35 + t158 * t24) + t156 * t23), t156 * (t158 * (-qJ(4) * t25 - t154 * t8 + t157 * t17) - t155 * (-pkin(3) * t25 - pkin(4) * t54 - pkin(6) * t38 - t230) - qJ(3) * t18) + t159 * (-pkin(2) * t18 + pkin(3) * t37 - qJ(4) * t26 - t154 * t17 - t157 * t8) - pkin(1) * t18 + qJ(2) * (t159 * (t155 * t37 + t158 * t26) + t156 * t25), t156 * (t158 * (-qJ(4) * t19 + t157 * t4 + t27 * t237) - t155 * (-pkin(3) * t19 - pkin(4) * t58 - pkin(6) * t29 - t6) - qJ(3) * t10) + t159 * (-pkin(2) * t10 - qJ(4) * t20 - t154 * t4 + t198 * t27) - pkin(1) * t10 + qJ(2) * (t159 * (t155 * t27 + t158 * t20) + t156 * t19), t156 * (t158 * (-qJ(4) * t2 + (-pkin(6) * t157 + t237) * t5) - t155 * (-pkin(3) * t2 + pkin(4) * t21 - pkin(6) * t6) - qJ(3) * t1) + t159 * (-pkin(2) * t1 - qJ(4) * t3 + (pkin(6) * t154 + t198) * t5) - pkin(1) * t1 + qJ(2) * (t159 * (t155 * t5 + t158 * t3) + t156 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, t209, -t133, t126, 0, 0, 0, 0, 0, 0, t94, t95, t91, t34, 0, 0, 0, 0, 0, 0, t40, t43, t39, t9, 0, 0, 0, 0, 0, 0, t14, t18, t10, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t119, -t129, t81, 0, 0, 0, 0, 0, 0, t56, t69, t59, t15, 0, 0, 0, 0, 0, 0, t23, t25, t19, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t89, t76, t45, 0, 0, 0, 0, 0, 0, t35, t37, t27, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, t53, -t75, -t49, t111, -t11, -t12, 0, 0;];
tauJ_reg = t32;
