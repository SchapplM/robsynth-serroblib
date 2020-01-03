% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:25
% EndTime: 2020-01-03 11:50:34
% DurationCPUTime: 2.95s
% Computational Cost: add. (7259->296), mult. (17799->377), div. (0->0), fcn. (12163->8), ass. (0->183)
t232 = 2 * qJD(2);
t152 = cos(pkin(8));
t195 = t152 * qJDD(1);
t138 = -qJDD(3) + t195;
t131 = -qJDD(4) + t138;
t153 = sin(qJ(4));
t154 = sin(qJ(3));
t156 = cos(qJ(4));
t157 = cos(qJ(3));
t151 = sin(pkin(8));
t203 = qJD(1) * t151;
t109 = (-t153 * t157 - t154 * t156) * t203;
t201 = qJD(1) * t157;
t181 = t151 * t201;
t182 = t154 * t203;
t111 = -t153 * t182 + t156 * t181;
t209 = t111 * t109;
t237 = t131 - t209;
t246 = t237 * pkin(4);
t245 = t153 * t237;
t244 = t156 * t237;
t202 = qJD(1) * t152;
t139 = -qJD(3) + t202;
t132 = -qJD(4) + t139;
t210 = t109 * t132;
t118 = (qJD(3) * t201 + qJDD(1) * t154) * t151;
t196 = t151 * qJDD(1);
t135 = t157 * t196;
t119 = -qJD(3) * t182 + t135;
t169 = -t118 * t153 + t119 * t156;
t68 = qJD(4) * t109 + t169;
t59 = t210 + t68;
t239 = qJ(5) * t59;
t145 = t151 ^ 2;
t159 = qJD(1) ^ 2;
t206 = t145 * t159;
t188 = t157 * t206;
t172 = -pkin(2) * t152 - pkin(6) * t151;
t125 = t172 * qJD(1);
t155 = sin(qJ(1));
t158 = cos(qJ(1));
t170 = g(2) * t155 - g(3) * t158;
t123 = -pkin(1) * t159 + qJDD(1) * qJ(2) - t170;
t238 = qJD(1) * t232 + t123;
t175 = -g(1) * t151 + t152 * t238;
t77 = t125 * t202 + t175;
t147 = t159 * qJ(2);
t150 = qJDD(1) * pkin(1);
t171 = -g(2) * t158 - g(3) * t155;
t120 = qJDD(2) - t147 - t150 - t171;
t99 = qJDD(1) * t172 + t120;
t89 = t157 * t99;
t43 = -pkin(3) * t138 - pkin(7) * t119 + t89 + (pkin(7) * t139 * t203 - pkin(3) * t188 - t77) * t154;
t115 = -pkin(3) * t139 - pkin(7) * t181;
t137 = t154 ^ 2 * t206;
t62 = t154 * t99 + t157 * t77;
t44 = -pkin(3) * t137 - pkin(7) * t118 + t115 * t139 + t62;
t24 = t153 * t44 - t156 * t43;
t240 = 2 * qJD(5);
t163 = t111 * t240 + t239 + t24 + t246;
t162 = -t163 - t246;
t107 = t109 ^ 2;
t130 = t132 ^ 2;
t71 = -t130 - t107;
t47 = t153 * t71 - t244;
t46 = pkin(3) * t47;
t243 = t162 + t46;
t178 = -t118 * t156 - t119 * t153;
t166 = qJD(4) * t111 - t178;
t25 = t153 * t43 + t156 * t44;
t87 = -pkin(4) * t132 - qJ(5) * t111;
t14 = -pkin(4) * t107 - qJ(5) * t166 + t109 * t240 + t132 * t87 + t25;
t108 = t111 ^ 2;
t85 = -t108 - t130;
t167 = pkin(4) * t85 - t14;
t72 = t131 + t209;
t215 = t153 * t72;
t52 = t156 * t85 + t215;
t51 = pkin(3) * t52;
t242 = t51 + t167;
t146 = t152 ^ 2;
t241 = t145 + t146;
t124 = t139 * t181;
t92 = t124 - t118;
t236 = t151 * t92;
t235 = t147 * t241 + t120 - t150;
t234 = t182 * (qJD(3) - t139) - t135;
t233 = pkin(4) * t166 - qJ(5) * t107 + t111 * t87 + qJDD(5);
t136 = t139 ^ 2;
t7 = t153 * t25 - t156 * t24;
t231 = pkin(3) * t7;
t48 = t156 * t71 + t245;
t30 = t154 * t48 + t157 * t47;
t230 = pkin(2) * t30;
t213 = t156 * t72;
t53 = -t153 * t85 + t213;
t32 = t154 * t53 + t157 * t52;
t229 = pkin(2) * t32;
t57 = (-qJD(4) - t132) * t111 + t178;
t35 = t153 * t57 - t156 * t59;
t36 = t153 * t59 + t156 * t57;
t16 = t154 * t36 + t157 * t35;
t228 = pkin(6) * t16;
t227 = pkin(6) * t30;
t226 = pkin(6) * t32;
t225 = pkin(7) * t35;
t224 = pkin(7) * t47;
t223 = pkin(7) * t52;
t222 = g(1) * t152;
t221 = t157 * t7;
t69 = -t107 - t108;
t220 = qJ(2) * (t152 * (-t154 * t35 + t157 * t36) + t151 * t69) - pkin(1) * t16;
t98 = t132 * t111;
t55 = t166 - t98;
t219 = qJ(2) * (t152 * (-t154 * t47 + t157 * t48) + t151 * t55) - pkin(1) * t30;
t60 = -(-qJD(4) + t132) * t109 + t169;
t218 = qJ(2) * (t152 * (-t154 * t52 + t157 * t53) + t151 * t60) - pkin(1) * t32;
t174 = pkin(3) * t118 - pkin(7) * t137 + t222;
t177 = -t115 * t157 - t125;
t49 = (t123 + (t232 - t177) * qJD(1)) * t151 + t174;
t216 = t153 * t49;
t214 = t156 * t49;
t212 = qJ(5) * t153;
t211 = qJ(5) * t156;
t208 = t132 * t153;
t207 = t132 * t156;
t129 = t154 * t188;
t116 = -t129 + t138;
t205 = t154 * t116;
t117 = -t129 - t138;
t204 = t157 * t117;
t194 = t51 - t25;
t11 = pkin(4) * t163;
t3 = t14 * t153 - t156 * t163;
t192 = pkin(3) * t3 - t11;
t189 = t157 ^ 2 * t206;
t187 = t152 * t209;
t34 = pkin(3) * t35;
t186 = -pkin(2) * t16 - t34;
t185 = -pkin(3) * t69 + pkin(7) * t36;
t184 = -pkin(3) * t55 + pkin(7) * t48;
t183 = -pkin(3) * t60 + pkin(7) * t53;
t180 = t151 * (t151 * t238 + t222) + t152 * t175;
t8 = t153 * t24 + t156 * t25;
t61 = t154 * t77 - t89;
t176 = t24 - t46;
t37 = t154 * t62 - t157 * t61;
t168 = -pkin(1) + t172;
t28 = t49 + t233;
t144 = t146 * qJDD(1);
t143 = t145 * qJDD(1);
t126 = t241 * t159;
t122 = -t137 + t189;
t121 = -t137 - t136;
t104 = -t189 - t136;
t94 = -t108 + t130;
t93 = t107 - t130;
t91 = t124 + t118;
t90 = -t135 + (qJD(3) + t139) * t182;
t84 = t121 * t154 + t204;
t79 = t108 - t107;
t78 = t104 * t157 + t205;
t76 = t222 + (t123 + (t232 + t125) * qJD(1)) * t151;
t63 = -t154 * t91 + t157 * t90;
t58 = t68 - t210;
t56 = t166 + t98;
t54 = pkin(4) * t59;
t39 = -pkin(4) * t60 + qJ(5) * t72;
t38 = t152 * t131 + t151 * (t157 * (-t109 * t156 - t111 * t153) - t154 * (-t109 * t153 + t111 * t156)) * t132;
t27 = t151 * (t157 * (t111 * t208 + t156 * t68) - t154 * (-t111 * t207 + t153 * t68)) + t187;
t26 = t151 * (t157 * (t109 * t207 + t153 * t166) - t154 * (t109 * t208 - t156 * t166)) - t187;
t22 = -qJ(5) * t85 + t28;
t21 = t151 * (t157 * (t156 * t93 + t215) - t154 * (t153 * t93 - t213)) + t152 * t56;
t20 = t151 * (t157 * (-t153 * t94 - t244) - t154 * (t156 * t94 - t245)) - t152 * t59;
t17 = -pkin(4) * t55 + qJ(5) * t71 - t174 + (qJD(1) * t177 - t238) * t151 - t233;
t13 = t151 * (t157 * (-t153 * t58 - t156 * t55) - t154 * (-t153 * t55 + t156 * t58)) - t152 * t79;
t9 = t163 + t239;
t6 = -pkin(4) * t69 + qJ(5) * t57 + t14;
t5 = -pkin(4) * t28 + qJ(5) * t14;
t4 = t14 * t156 + t153 * t163;
t2 = t154 * t8 + t221;
t1 = t154 * t4 + t157 * t3;
t10 = [0, 0, 0, 0, 0, qJDD(1), t171, t170, 0, 0, t143, 0.2e1 * t151 * t195, 0, t144, 0, 0, -t235 * t152, t235 * t151, pkin(1) * t126 + qJ(2) * (t144 + t143) + t180, -pkin(1) * t120 + qJ(2) * t180, (t151 * (t139 * t182 + t119) - t152 * t154 * t206) * t157, t151 * (t154 * t234 + t157 * t92) - t152 * t122, t151 * (t204 - t154 * (t136 - t189)) + t152 * t90, (t152 * t188 - t236) * t154, t151 * (t157 * (t137 - t136) + t205) + t152 * t91, t152 * t138, t151 * (-pkin(6) * t84 + t154 * t76) + t152 * (-pkin(2) * t84 + t61) - pkin(1) * t84 + qJ(2) * (t152 * (-t117 * t154 + t121 * t157) - t236), t151 * (-pkin(6) * t78 + t157 * t76) + t152 * (-pkin(2) * t78 + t62) - pkin(1) * t78 + qJ(2) * (t152 * (-t104 * t154 + t116 * t157) - t234 * t151), -t151 * t37 + qJ(2) * (t152 * (-t154 * t90 - t157 * t91) - t151 * (t137 + t189)) + t168 * t63, qJ(2) * (t152 * (t154 * t61 + t157 * t62) + t151 * t76) + t168 * t37, t27, t13, t20, t26, t21, t38, t151 * (t157 * (t216 - t224) - t154 * (t184 - t214) - t227) + t152 * (t176 - t230) + t219, t151 * (t157 * (t214 - t223) - t154 * (t183 + t216) - t226) + t152 * (-t194 - t229) + t218, t151 * (t157 * (-t7 - t225) - t154 * (t185 + t8) - t228) + t152 * t186 + t220, t151 * (-pkin(7) * t221 - t154 * (-pkin(3) * t49 + pkin(7) * t8) - pkin(6) * t2) + t152 * (-pkin(2) * t2 - t231) - pkin(1) * t2 + qJ(2) * (t152 * (-t154 * t7 + t157 * t8) + t151 * t49), t27, t13, t20, t26, t21, t38, t151 * (t157 * (-t153 * t17 + t211 * t237 - t224) - t154 * (t156 * t17 + t212 * t237 + t184) - t227) + t152 * (-t230 - t243) + t219, t151 * (t157 * (-t153 * t39 + t156 * t22 - t223) - t154 * (t153 * t22 + t156 * t39 + t183) - t226) + t152 * (-t229 - t242) + t218, t151 * (t157 * (-t153 * t6 + t156 * t9 - t225) - t154 * (t153 * t9 + t156 * t6 + t185) - t228) + t152 * (t186 + t54) + t220, t151 * (t157 * (-pkin(7) * t3 - t153 * t5 + t163 * t211) - t154 * (-pkin(3) * t28 + pkin(7) * t4 + t156 * t5 + t163 * t212) - pkin(6) * t1) + t152 * (-pkin(2) * t1 - t192) - pkin(1) * t1 + qJ(2) * (t152 * (-t154 * t3 + t157 * t4) + t151 * t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t195, t196, -t126, t120, 0, 0, 0, 0, 0, 0, t84, t78, t63, t37, 0, 0, 0, 0, 0, 0, t30, t32, t16, t2, 0, 0, 0, 0, 0, 0, t30, t32, t16, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t122, -t90, -t129, -t91, -t138, -t61, -t62, 0, 0, -t209, t79, t59, t209, -t56, -t131, -t176, t194, t34, t231, -t209, t79, t59, t209, -t56, -t131, t243, t242, -t54 + t34, t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t79, t59, t209, -t56, -t131, -t24, -t25, 0, 0, -t209, t79, t59, t209, -t56, -t131, t162, t167, -t54, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t60, t69, t28;];
tauJ_reg = t10;
