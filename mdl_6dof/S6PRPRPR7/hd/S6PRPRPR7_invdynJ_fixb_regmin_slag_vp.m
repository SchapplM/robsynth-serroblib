% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% tau_reg [6x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:50
% EndTime: 2019-03-08 19:53:56
% DurationCPUTime: 2.72s
% Computational Cost: add. (1531->348), mult. (3183->452), div. (0->0), fcn. (2360->10), ass. (0->193)
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t109 = cos(pkin(6));
t208 = qJD(1) * t109;
t117 = -pkin(2) - pkin(8);
t115 = cos(qJ(2));
t107 = sin(pkin(6));
t209 = qJD(1) * t107;
t175 = t115 * t209;
t155 = qJD(3) - t175;
t55 = t117 * qJD(2) + t155;
t29 = t111 * t208 - t114 * t55;
t257 = qJD(5) + t29;
t173 = qJD(4) * t208;
t193 = qJDD(1) * t109;
t201 = qJD(4) * t114;
t194 = qJDD(1) * t107;
t171 = t115 * t194;
t112 = sin(qJ(2));
t176 = t112 * t209;
t74 = qJD(2) * t176;
t141 = qJDD(3) + t74 - t171;
t31 = t117 * qJDD(2) + t141;
t168 = -t114 * t193 - t55 * t201 + (t173 - t31) * t111;
t191 = qJDD(4) * qJ(5);
t128 = qJD(4) * qJD(5) - t168 + t191;
t16 = -qJD(4) * pkin(4) + t257;
t203 = qJD(4) * qJ(5);
t30 = t111 * t55 + t114 * t208;
t17 = -t203 - t30;
t202 = qJD(4) * t111;
t186 = t111 * t193 + t114 * t173 + t55 * t202;
t248 = -t114 * t31 + qJDD(5);
t138 = t186 + t248;
t227 = qJDD(4) * pkin(4);
t7 = t138 - t227;
t124 = (t111 * t16 - t114 * t17) * qJD(4) + t128 * t111 - t7 * t114;
t106 = sin(pkin(10));
t108 = cos(pkin(10));
t221 = t109 * t115;
t49 = t106 * t112 - t108 * t221;
t51 = t106 * t221 + t108 * t112;
t162 = -g(1) * t51 - g(2) * t49;
t223 = t107 * t115;
t150 = -g(3) * t223 - t162;
t256 = t124 - t150;
t195 = qJD(2) * qJD(4);
t172 = t114 * t195;
t187 = t111 * qJDD(2);
t255 = t172 + t187;
t110 = sin(qJ(6));
t199 = t110 * qJD(4);
t113 = cos(qJ(6));
t205 = qJD(2) * t113;
t254 = t111 * t205 - t199;
t119 = qJD(2) ^ 2;
t216 = t115 * t119;
t146 = qJDD(2) * t112 + t216;
t225 = t107 * t112;
t178 = qJD(2) * t225;
t54 = t109 * t114 - t111 * t223;
t23 = t54 * qJD(4) - t114 * t178;
t53 = t109 * t111 + t114 * t223;
t253 = -t23 * qJD(4) - t53 * qJDD(4) + (t111 * t146 + t112 * t172) * t107;
t102 = t114 * qJDD(2);
t170 = t111 * t195;
t136 = t170 - t102;
t22 = t53 * qJD(4) - t111 * t178;
t252 = -t22 * qJD(4) + t54 * qJDD(4) + (t112 * t136 - t114 * t216) * t107;
t251 = pkin(4) * t201 + qJ(5) * t202;
t197 = t114 * qJD(2);
t94 = qJD(6) + t197;
t234 = t110 * t94;
t61 = -qJDD(6) + t136;
t46 = t113 * t61;
t142 = -qJD(6) * t234 - t46;
t222 = t109 * t112;
t50 = t106 * t115 + t108 * t222;
t52 = -t106 * t222 + t108 * t115;
t250 = -g(1) * t52 - g(2) * t50;
t249 = t94 - qJD(6);
t213 = pkin(5) * t197 + t257;
t116 = -pkin(4) - pkin(9);
t206 = qJD(2) * t111;
t15 = -pkin(5) * t206 + t30;
t13 = t15 + t203;
t247 = -t116 * t61 + (t13 - t15) * t94;
t198 = t113 * qJD(4);
t65 = t110 * t206 + t198;
t12 = t65 * qJD(6) + t110 * qJDD(4) - t113 * t255;
t188 = qJDD(4) * t117;
t207 = qJD(2) * qJ(3);
t67 = t176 + t207;
t244 = (t176 - t67 - t207) * qJD(4) - t188;
t160 = pkin(4) * t206 + t176;
t217 = t114 * qJ(5);
t165 = qJ(3) - t217;
t34 = t165 * qJD(2) + t160;
t103 = t111 * pkin(4);
t68 = t103 + t165;
t243 = (-qJD(2) * t68 + t176 - t34) * qJD(4) - t188;
t3 = -pkin(5) * t255 + t128;
t240 = t3 * t111;
t239 = t254 * t94;
t238 = t65 * t94;
t237 = pkin(5) - t117;
t11 = qJD(6) * t254 + t113 * qJDD(4) + t110 * t255;
t236 = t11 * t113;
t235 = t110 * t61;
t233 = t111 * t13;
t232 = t111 * t65;
t66 = pkin(4) * t197 + qJ(5) * t206;
t230 = qJD(4) * t254;
t229 = qJD(4) * t65;
t228 = qJDD(2) * pkin(2);
t226 = t107 * t111;
t224 = t107 * t114;
t220 = t110 * t114;
t219 = t111 * t112;
t218 = t113 * t114;
t215 = t30 * qJD(4);
t214 = t67 * qJD(2);
t104 = t111 ^ 2;
t105 = t114 ^ 2;
t212 = t104 - t105;
t211 = t104 + t105;
t118 = qJD(4) ^ 2;
t210 = t118 + t119;
t204 = qJD(2) * t115;
t200 = qJD(6) * t113;
t196 = t114 * qJD(5);
t192 = qJDD(2) * qJ(3);
t190 = qJDD(4) * t111;
t189 = qJDD(4) * t114;
t185 = g(3) * (pkin(2) * t223 + qJ(3) * t225);
t184 = t94 * t199;
t183 = t94 * t198;
t180 = t114 * t119 * t111;
t2 = -t136 * pkin(5) + t116 * qJDD(4) + t138;
t135 = qJD(3) + (qJD(4) * pkin(9) - qJD(5)) * t114;
t144 = t111 * pkin(9) + t165;
t82 = t112 * t194;
t167 = pkin(4) * t255 + qJ(5) * t170 + t82;
t8 = t144 * qJDD(2) + (t135 + t175) * qJD(2) + t167;
t179 = -t110 * t8 + t113 * t2;
t177 = t107 * t204;
t169 = qJD(4) * t237;
t164 = t211 * qJDD(2);
t156 = -t217 + t103;
t154 = qJD(3) + t175;
t10 = t116 * qJD(4) + t213;
t24 = t144 * qJD(2) + t160;
t4 = t113 * t10 - t110 * t24;
t5 = t110 * t10 + t113 * t24;
t152 = t112 * (-qJD(2) * pkin(2) + t155) + t115 * t67;
t143 = t94 * t200 - t235;
t26 = t53 * t110 + t113 * t225;
t25 = -t110 * t225 + t53 * t113;
t140 = t110 * t115 + t112 * t218;
t139 = -t112 * t220 + t113 * t115;
t19 = t106 * t224 + t51 * t111;
t21 = t108 * t224 - t49 * t111;
t137 = -g(1) * t19 + g(2) * t21 - g(3) * t54;
t132 = g(3) * t225 - t250;
t131 = t132 - t82;
t130 = t150 + t171;
t18 = t106 * t226 - t51 * t114;
t20 = t108 * t226 + t49 * t114;
t129 = g(1) * t18 - g(2) * t20 + g(3) * t53 - t186;
t127 = t137 - t168;
t126 = -t117 * t118 - t132;
t125 = qJDD(3) - t130;
t123 = t3 + (pkin(9) * t197 - qJD(6) * t116 + t66) * t94 + t137;
t122 = t34 * t197 - t129 + t248;
t42 = qJD(3) - t196 + t251;
t9 = t165 * qJDD(2) + (t154 - t196) * qJD(2) + t167;
t121 = -qJDD(2) * t68 - t9 + (-t42 + t175) * qJD(2) - t126;
t32 = t154 * qJD(2) + t192 + t82;
t120 = t155 * qJD(2) + t126 + t192 + t32;
t70 = t237 * t114;
t69 = t237 * t111;
t60 = t114 * t169;
t59 = t111 * t169;
t58 = t103 + t144;
t57 = -t210 * t111 + t189;
t56 = t210 * t114 + t190;
t48 = t146 * t107;
t47 = (-qJDD(2) * t115 + t112 * t119) * t107;
t45 = t51 * pkin(2);
t44 = t49 * pkin(2);
t35 = t141 - t228;
t33 = t135 + t251;
t1 = [qJDD(1) - g(3), 0, -t47, -t48, t47, t48, t109 ^ 2 * qJDD(1) - g(3) + (t152 * qJD(2) + t112 * t32 - t115 * t35) * t107, 0, 0, 0, 0, 0, t253, -t252 (-t111 * t54 + t114 * t53) * qJDD(2) + (t111 * t22 + t114 * t23 + (-t111 * t53 - t114 * t54) * qJD(4)) * qJD(2), -t253, t252, t16 * t23 + t17 * t22 + t7 * t53 + t128 * t54 - g(3) + (t112 * t9 + t204 * t34) * t107, 0, 0, 0, 0, 0 (-qJD(6) * t26 - t110 * t177 + t23 * t113) * t94 - t25 * t61 + t22 * t254 + t54 * t12 -(qJD(6) * t25 + t23 * t110 + t113 * t177) * t94 + t26 * t61 - t22 * t65 + t54 * t11; 0, qJDD(2), t130, t131, t125 - 0.2e1 * t228, 0.2e1 * qJD(2) * qJD(3) - t131 + 0.2e1 * t192, t32 * qJ(3) + t67 * qJD(3) - t35 * pkin(2) - g(1) * (t52 * qJ(3) - t45) - g(2) * (t50 * qJ(3) - t44) - t185 - t152 * t209, t105 * qJDD(2) - 0.2e1 * t114 * t170, -0.2e1 * t111 * t102 + 0.2e1 * t212 * t195, -t118 * t111 + t189, -t118 * t114 - t190, 0, t120 * t111 - t244 * t114, t244 * t111 + t120 * t114, -t117 * t164 + t211 * t74 - t256, t121 * t111 + t243 * t114, -t243 * t111 + t121 * t114, t9 * t68 + t34 * t42 - g(1) * (-t51 * pkin(8) - t45) - g(2) * (-t49 * pkin(8) - t44) - t185 + t124 * t117 + (-g(3) * (pkin(8) * t115 + t112 * t156) + (-t34 * t115 + (t111 * t17 + t114 * t16) * t112) * qJD(1)) * t107 + t250 * (qJ(3) + t156) t200 * t232 + (t111 * t11 + t201 * t65) * t110 (t110 * t254 + t113 * t65) * t201 + (t236 - t110 * t12 + (-t110 * t65 + t113 * t254) * qJD(6)) * t111 (t11 + t184) * t114 + (t143 - t229) * t111 (-t12 + t183) * t114 + (t142 - t230) * t111, -t61 * t114 - t202 * t94 (-t110 * t33 - t113 * t59) * t94 - (-t110 * t58 + t113 * t70) * t61 + t179 * t114 + t60 * t254 - t69 * t12 - t113 * t240 - g(1) * (-t51 * t113 - t220 * t52) - g(2) * (-t49 * t113 - t220 * t50) + (-t4 * t111 - t13 * t218) * qJD(4) + ((-t110 * t70 - t113 * t58) * t94 - t5 * t114 + t110 * t233) * qJD(6) + (-g(3) * t139 + (t140 * t94 + t219 * t254) * qJD(1)) * t107, t5 * t202 - t69 * t11 - t60 * t65 + (-(qJD(6) * t70 + t33) * t94 + t58 * t61 + qJD(6) * t233 + (-qJD(6) * t10 - t250 - t8) * t114) * t113 + (-(-qJD(6) * t58 - t59) * t94 + t70 * t61 + t240 + (t13 * qJD(4) + qJD(6) * t24 - t2) * t114 + t162) * t110 + (g(3) * t140 + (t139 * t94 - t219 * t65) * qJD(1)) * t107; 0, 0, 0, 0, qJDD(2), -t119, t125 - t214 + t74 - t228, 0, 0, 0, 0, 0, t57, -t56, -t164, -t57, t56, -t34 * qJD(2) + t256, 0, 0, 0, 0, 0, qJD(2) * t234 + (t12 + t183) * t111 + (-t142 - t230) * t114, t94 * t205 + (t11 - t184) * t111 + (t143 + t229) * t114; 0, 0, 0, 0, 0, 0, 0, t180, -t212 * t119, t102, -t187, qJDD(4), t215 + (t31 - t214) * t114 + t129, -t29 * qJD(4) + t206 * t67 - t127 (-pkin(4) * t114 - qJ(5) * t111) * qJDD(2), t206 * t66 + t122 - t215 - 0.2e1 * t227, 0.2e1 * t191 + (0.2e1 * qJD(5) + t29) * qJD(4) + (-t111 * t34 + t114 * t66) * qJD(2) + t127, t128 * qJ(5) - t7 * pkin(4) - t34 * t66 - t16 * t30 - g(1) * (-t18 * pkin(4) + t19 * qJ(5)) - g(2) * (t20 * pkin(4) - t21 * qJ(5)) - g(3) * (-t53 * pkin(4) + t54 * qJ(5)) - t257 * t17, -t234 * t65 + t236 (-t12 - t238) * t113 + (-t11 - t239) * t110 (-t220 * t94 + t232) * qJD(2) + t142 (t111 * t254 - t218 * t94) * qJD(2) - t143, t94 * t206, qJ(5) * t12 + t123 * t110 + t247 * t113 + t4 * t206 - t213 * t254, qJ(5) * t11 - t247 * t110 + t123 * t113 - t5 * t206 + t213 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, qJDD(4) - t180, -t105 * t119 - t118, t17 * qJD(4) + t122 - t227, 0, 0, 0, 0, 0, -t234 * t94 + t230 - t46, -t113 * t94 ^ 2 - t229 + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 * t254, -t254 ^ 2 + t65 ^ 2, t11 - t239, -t12 + t238, -t61, -t13 * t65 - g(1) * (-t52 * t110 + t18 * t113) - g(2) * (-t50 * t110 - t20 * t113) - g(3) * t25 + t179 + t249 * t5, -t113 * t8 - t110 * t2 - t13 * t254 - g(1) * (-t18 * t110 - t52 * t113) - g(2) * (t20 * t110 - t50 * t113) + g(3) * t26 + t249 * t4;];
tau_reg  = t1;
