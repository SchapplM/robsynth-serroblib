% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRR11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:59
% EndTime: 2024-09-27 21:46:02
% DurationCPUTime: 1.31s
% Computational Cost: add. (14372->343), mult. (36315->513), div. (0->0), fcn. (28288->12), ass. (0->210)
t172 = sin(qJ(5));
t171 = cos(pkin(5));
t161 = t171 * qJDD(2) + qJDD(3);
t157 = qJDD(4) + t161;
t148 = qJDD(5) + t157;
t173 = sin(qJ(4));
t177 = cos(qJ(4));
t169 = sin(pkin(5));
t178 = cos(qJ(3));
t213 = qJD(2) * t178;
t205 = t169 * t213;
t174 = sin(qJ(3));
t214 = qJD(2) * t169;
t206 = t174 * t214;
t132 = t173 * t206 - t177 * t205;
t134 = (t178 * t173 + t174 * t177) * t214;
t176 = cos(qJ(5));
t109 = t176 * t132 + t134 * t172;
t111 = -t132 * t172 + t134 * t176;
t80 = t111 * t109;
t243 = -t80 + t148;
t250 = t172 * t243;
t114 = t134 * t132;
t240 = -t114 + t157;
t249 = t173 * t240;
t248 = t176 * t243;
t247 = t177 * t240;
t216 = -g(3) + qJDD(1);
t180 = qJD(2) ^ 2;
t168 = sin(pkin(10));
t170 = cos(pkin(10));
t153 = g(1) * t168 - g(2) * t170;
t154 = -g(1) * t170 - g(2) * t168;
t175 = sin(qJ(2));
t179 = cos(qJ(2));
t201 = t179 * t153 - t154 * t175;
t237 = t169 * pkin(7);
t117 = qJDD(2) * pkin(2) + t180 * t237 + t201;
t220 = t171 * t117;
t246 = t169 * t216 + t220;
t162 = qJD(2) * t171 + qJD(3);
t159 = qJD(4) + t162;
t152 = qJD(5) + t159;
t100 = t152 * t109;
t211 = qJDD(2) * t174;
t140 = (qJD(3) * t213 + t211) * t169;
t212 = qJDD(2) * t169;
t141 = -qJD(3) * t206 + t178 * t212;
t202 = t140 * t173 - t177 * t141;
t95 = -qJD(4) * t134 - t202;
t96 = -t132 * qJD(4) + t177 * t140 + t173 * t141;
t59 = -t109 * qJD(5) + t172 * t95 + t176 * t96;
t245 = -t100 + t59;
t125 = t159 * t132;
t244 = -t125 + t96;
t242 = -t96 - t125;
t165 = t169 ^ 2;
t215 = qJD(2) * t162;
t241 = t165 * (-t171 * t180 + t215);
t146 = t162 * t205;
t120 = -t146 + t140;
t160 = t162 ^ 2;
t167 = t178 ^ 2;
t221 = t165 * t180;
t209 = t167 * t221;
t142 = -t160 - t209;
t150 = t174 * t178 * t221;
t138 = t161 + t150;
t107 = t109 ^ 2;
t108 = t111 ^ 2;
t130 = t132 ^ 2;
t131 = t134 ^ 2;
t149 = t152 ^ 2;
t156 = t159 ^ 2;
t185 = -t153 * t175 - t154 * t179;
t118 = -pkin(2) * t180 + pkin(7) * t212 - t185;
t89 = t118 * t174 - t246 * t178;
t68 = t138 * pkin(3) - t120 * pkin(8) - t89;
t218 = t178 * t118;
t69 = t141 * pkin(8) + t218 + t142 * pkin(3) + (t220 + (pkin(8) * t215 + t216) * t169) * t174;
t39 = t173 * t69 - t177 * t68;
t40 = t173 * t68 + t177 * t69;
t20 = t173 * t40 - t177 * t39;
t239 = pkin(3) * t20;
t181 = (-qJD(4) + t159) * t134 - t202;
t55 = t173 * t181 + t177 * t242;
t238 = pkin(3) * t55;
t31 = t240 * pkin(4) + t242 * pkin(9) - t39;
t199 = pkin(4) * t159 - pkin(9) * t134;
t33 = -t130 * pkin(4) + t95 * pkin(9) - t159 * t199 + t40;
t16 = t172 * t33 - t176 * t31;
t17 = t172 * t31 + t176 * t33;
t8 = -t16 * t176 + t17 * t172;
t236 = t173 * t8;
t235 = t177 * t8;
t112 = t117 * t169 - t171 * t216;
t79 = pkin(3) * t141 + pkin(8) * t209 + t112 - (pkin(3) * t162 - pkin(8) * t206) * t206;
t51 = pkin(4) * t95 + pkin(9) * t130 - t134 * t199 + t79;
t234 = t172 * t51;
t76 = t80 + t148;
t233 = t172 * t76;
t232 = t173 * t79;
t231 = t174 * t20;
t230 = t176 * t51;
t229 = t176 * t76;
t228 = t177 * t79;
t104 = t114 + t157;
t227 = t104 * t173;
t226 = t104 * t177;
t225 = t152 * t172;
t224 = t152 * t176;
t223 = t159 * t173;
t222 = t159 * t177;
t219 = t174 * t138;
t139 = -t150 + t161;
t217 = t178 * t139;
t166 = t174 ^ 2;
t210 = t166 * t221;
t208 = t171 * t80;
t207 = t171 * t114;
t9 = t16 * t172 + t176 * t17;
t204 = t172 * t96 - t176 * t95;
t21 = t173 * t39 + t177 * t40;
t74 = -t149 - t107;
t53 = t172 * t74 + t248;
t200 = pkin(4) * t53 - t16;
t2 = t173 * t9 + t235;
t3 = t177 * t9 - t236;
t198 = t174 * t3 + t178 * t2;
t182 = (-qJD(5) + t152) * t111 - t204;
t49 = t100 + t59;
t23 = t172 * t182 - t176 * t49;
t25 = t172 * t49 + t176 * t182;
t13 = t173 * t25 + t177 * t23;
t14 = -t173 * t23 + t177 * t25;
t197 = t13 * t178 + t14 * t174;
t196 = t174 * t21 + t178 * t20;
t54 = t176 * t74 - t250;
t28 = t173 * t54 + t177 * t53;
t29 = -t173 * t53 + t177 * t54;
t195 = t174 * t29 + t178 * t28;
t91 = -t108 - t149;
t60 = t176 * t91 - t233;
t61 = -t172 * t91 - t229;
t36 = t173 * t61 + t177 * t60;
t37 = -t173 * t60 + t177 * t61;
t194 = t174 * t37 + t178 * t36;
t56 = -t173 * t242 + t177 * t181;
t193 = t174 * t56 + t178 * t55;
t102 = -t156 - t130;
t72 = t102 * t173 + t247;
t73 = t102 * t177 - t249;
t192 = t174 * t73 + t178 * t72;
t116 = -t131 - t156;
t81 = t116 * t177 - t227;
t82 = -t116 * t173 - t226;
t191 = t174 * t82 + t178 * t81;
t90 = t246 * t174 + t218;
t190 = t174 * t90 - t178 * t89;
t189 = t174 * t89 + t178 * t90;
t145 = t162 * t206;
t121 = t141 + t145;
t188 = -t120 * t178 + t121 * t174;
t128 = -t210 - t160;
t187 = t128 * t178 - t139 * t174;
t186 = t138 * t178 + t142 * t174;
t183 = pkin(4) * t60 - t17;
t144 = (-t166 - t167) * t221;
t143 = (t166 - t167) * t221;
t124 = -t131 + t156;
t123 = t130 - t156;
t122 = t141 - t145;
t119 = (t211 + (qJD(3) + t162) * t213) * t169;
t113 = t131 - t130;
t99 = -t108 + t149;
t98 = t107 - t149;
t97 = -t130 - t131;
t83 = (qJD(4) + t159) * t134 + t202;
t78 = t108 - t107;
t71 = (-t109 * t176 + t111 * t172) * t152;
t70 = (-t109 * t172 - t111 * t176) * t152;
t66 = -t107 - t108;
t65 = t176 * t98 - t233;
t64 = -t172 * t99 + t248;
t63 = t172 * t98 + t229;
t62 = t176 * t99 + t250;
t58 = -qJD(5) * t111 - t204;
t45 = (qJD(5) + t152) * t111 + t204;
t44 = -t111 * t225 + t176 * t59;
t43 = t111 * t224 + t172 * t59;
t42 = t109 * t224 - t172 * t58;
t41 = t109 * t225 + t176 * t58;
t35 = pkin(3) * t81 - t40;
t34 = pkin(3) * t72 - t39;
t32 = -pkin(9) * t60 - t230;
t27 = -pkin(9) * t53 - t234;
t26 = -t172 * t245 - t176 * t45;
t24 = -t172 * t45 + t176 * t245;
t22 = pkin(4) * t23;
t19 = -pkin(4) * t245 + pkin(9) * t61 - t234;
t18 = -pkin(4) * t45 + pkin(9) * t54 + t230;
t12 = pkin(3) * t36 + t183;
t11 = pkin(3) * t13 + t22;
t10 = pkin(3) * t28 + t200;
t7 = pkin(4) * t8;
t6 = pkin(4) * t51 + pkin(9) * t9;
t5 = -pkin(9) * t23 - t8;
t4 = -pkin(4) * t66 + pkin(9) * t25 + t9;
t1 = pkin(3) * t2 + t7;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t216, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, 0, 0, 0, 0, 0, 0, -t122 * t171 + t169 * t186, t119 * t171 + t169 * t187, t144 * t171 + t169 * t188, -t112 * t171 + t169 * t190, 0, 0, 0, 0, 0, 0, t169 * t192 + t171 * t83, t169 * t191 + t171 * t244, t169 * t193 + t171 * t97, t169 * t196 - t171 * t79, 0, 0, 0, 0, 0, 0, t169 * t195 + t171 * t45, t169 * t194 + t171 * t245, t169 * t197 + t171 * t66, t169 * t198 - t171 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t201, t185, 0, 0, (t140 * t169 + t178 * t241) * t174, t171 * t143 + (t174 * t122 + (t140 + t146) * t178) * t169, t171 * t120 + (t219 + t178 * (t160 - t210)) * t169, (t141 * t169 - t174 * t241) * t178, t171 * t121 + (t174 * (-t160 + t209) + t217) * t169, t171 * t161, (pkin(2) * t186 - t89) * t171 + (t178 * t112 + pkin(2) * t122 + pkin(7) * (t142 * t178 - t219)) * t169, (pkin(2) * t187 - t90) * t171 + (-t174 * t112 - pkin(2) * t119 + pkin(7) * (-t128 * t174 - t217)) * t169, pkin(2) * t188 * t171 + (-pkin(2) * t144 + pkin(7) * (t120 * t174 + t121 * t178) + t189) * t169, pkin(2) * (t112 * t169 + t171 * t190) + t189 * t237, t207 + (t174 * (-t134 * t223 + t177 * t96) + t178 * (t134 * t222 + t173 * t96)) * t169, t171 * t113 + (t174 * (-t173 * t244 - t177 * t83) + t178 * (-t173 * t83 + t177 * t244)) * t169, -t171 * t242 + (t174 * (-t124 * t173 + t247) + t178 * (t124 * t177 + t249)) * t169, -t207 + (t174 * (t132 * t222 - t173 * t95) + t178 * (t132 * t223 + t177 * t95)) * t169, t171 * t181 + (t174 * (t123 * t177 - t227) + t178 * (t123 * t173 + t226)) * t169, t171 * t157 + (t174 * (-t132 * t177 + t134 * t173) + t178 * (-t132 * t173 - t134 * t177)) * t169 * t159, (pkin(2) * t192 + t34) * t171 + (t174 * (-pkin(8) * t72 - t232) + t178 * (-pkin(3) * t83 + pkin(8) * t73 + t228) - pkin(2) * t83 + pkin(7) * (-t174 * t72 + t178 * t73)) * t169, (pkin(2) * t191 + t35) * t171 + (t174 * (-pkin(8) * t81 - t228) + t178 * (-pkin(3) * t244 + pkin(8) * t82 - t232) - pkin(2) * t244 + pkin(7) * (-t174 * t81 + t178 * t82)) * t169, (pkin(2) * t193 + t238) * t171 + (t174 * (-pkin(8) * t55 - t20) + t178 * (-pkin(3) * t97 + pkin(8) * t56 + t21) - pkin(2) * t97 + pkin(7) * (-t174 * t55 + t178 * t56)) * t169, (pkin(2) * t196 + t239) * t171 + (-pkin(8) * t231 + t178 * (pkin(3) * t79 + pkin(8) * t21) + pkin(2) * t79 + pkin(7) * (t178 * t21 - t231)) * t169, t208 + (t174 * (-t173 * t43 + t177 * t44) + t178 * (t173 * t44 + t177 * t43)) * t169, t171 * t78 + (t174 * (-t173 * t24 + t177 * t26) + t178 * (t173 * t26 + t177 * t24)) * t169, t171 * t49 + (t174 * (-t173 * t62 + t177 * t64) + t178 * (t173 * t64 + t177 * t62)) * t169, -t208 + (t174 * (-t173 * t41 + t177 * t42) + t178 * (t173 * t42 + t177 * t41)) * t169, t171 * t182 + (t174 * (-t173 * t63 + t177 * t65) + t178 * (t173 * t65 + t177 * t63)) * t169, t171 * t148 + (t174 * (-t173 * t70 + t177 * t71) + t178 * (t173 * t71 + t177 * t70)) * t169, (pkin(2) * t195 + t10) * t171 + (t174 * (-pkin(8) * t28 - t173 * t18 + t177 * t27) + t178 * (-pkin(3) * t45 + pkin(8) * t29 + t173 * t27 + t177 * t18) - pkin(2) * t45 + pkin(7) * (-t174 * t28 + t178 * t29)) * t169, (pkin(2) * t194 + t12) * t171 + (t174 * (-pkin(8) * t36 - t173 * t19 + t177 * t32) + t178 * (-pkin(3) * t245 + pkin(8) * t37 + t173 * t32 + t177 * t19) - pkin(2) * t245 + pkin(7) * (-t174 * t36 + t178 * t37)) * t169, (pkin(2) * t197 + t11) * t171 + (t174 * (-pkin(8) * t13 - t173 * t4 + t177 * t5) + t178 * (-pkin(3) * t66 + pkin(8) * t14 + t173 * t5 + t177 * t4) - pkin(2) * t66 + pkin(7) * (-t13 * t174 + t14 * t178)) * t169, (pkin(2) * t198 + t1) * t171 + (t174 * (-pkin(8) * t2 - pkin(9) * t235 - t173 * t6) + t178 * (pkin(3) * t51 + pkin(8) * t3 - pkin(9) * t236 + t177 * t6) + pkin(2) * t51 + pkin(7) * (-t174 * t2 + t178 * t3)) * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, t143, t120, t150, t121, t161, -t89, -t90, 0, 0, t114, t113, -t242, -t114, t181, t157, t34, t35, t238, t239, t80, t78, t49, -t80, t182, t148, t10, t12, t11, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t113, -t242, -t114, t181, t157, -t39, -t40, 0, 0, t80, t78, t49, -t80, t182, t148, t200, t183, t22, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t78, t49, -t80, t182, t148, -t16, -t17, 0, 0;];
tauJ_reg = t15;
