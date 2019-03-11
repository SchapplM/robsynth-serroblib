% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% tau_reg [6x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:12
% EndTime: 2019-03-09 01:47:17
% DurationCPUTime: 2.09s
% Computational Cost: add. (2633->313), mult. (5016->400), div. (0->0), fcn. (3406->12), ass. (0->176)
t134 = sin(qJ(6));
t129 = sin(pkin(10));
t135 = sin(qJ(4));
t202 = qJD(1) * t135;
t137 = cos(qJ(4));
t215 = cos(pkin(10));
t179 = t215 * t137;
t89 = qJD(1) * t179;
t67 = t129 * t202 - t89;
t210 = qJD(6) - t67;
t136 = cos(qJ(6));
t236 = t136 * t210;
t190 = t135 * qJDD(1);
t180 = t215 * t135;
t73 = t129 * t137 + t180;
t68 = t73 * qJD(4);
t44 = qJD(1) * t68 - qJDD(1) * t179 + t129 * t190;
t43 = -qJDD(6) + t44;
t238 = t134 * t43 - t210 * t236;
t198 = t136 * qJD(4);
t69 = t73 * qJD(1);
t51 = -t134 * t69 - t198;
t237 = t210 * t51;
t200 = qJD(4) * t135;
t70 = qJD(4) * t179 - t129 * t200;
t160 = -t210 * t70 + t43 * t73;
t212 = qJD(6) * t73;
t235 = -t160 * t134 + t212 * t236;
t115 = t137 * qJDD(3);
t130 = sin(pkin(9));
t131 = cos(pkin(9));
t191 = t131 * qJDD(1);
t138 = -pkin(1) - pkin(2);
t85 = t138 * qJDD(1) + qJDD(2);
t223 = qJ(2) * t191 + t130 * t85;
t196 = qJD(1) * qJD(2);
t93 = t131 * t196;
t50 = t93 + t223;
t47 = -qJDD(1) * pkin(7) + t50;
t145 = qJ(5) * qJDD(1) + qJD(1) * qJD(5) - qJD(4) * qJD(3) - t47;
t203 = qJD(1) * t131;
t86 = t138 * qJD(1) + qJD(2);
t65 = qJ(2) * t203 + t130 * t86;
t60 = -qJD(1) * pkin(7) + t65;
t174 = qJ(5) * qJD(1) - t60;
t157 = t174 * qJD(4);
t11 = qJDD(4) * pkin(4) + t145 * t135 + t137 * t157 + t115;
t12 = (qJDD(3) + t157) * t135 - t145 * t137;
t3 = t215 * t11 - t129 * t12;
t1 = -qJDD(4) * pkin(5) - t3;
t42 = t135 * qJD(3) - t174 * t137;
t220 = t129 * t42;
t41 = qJ(5) * t202 + t137 * qJD(3) - t135 * t60;
t35 = qJD(4) * pkin(4) + t41;
t15 = t215 * t35 - t220;
t13 = -qJD(4) * pkin(5) - t15;
t211 = t129 * t135;
t152 = t179 - t211;
t120 = t137 * pkin(4);
t204 = qJD(1) * t130;
t64 = -qJ(2) * t204 + t131 * t86;
t59 = qJD(1) * pkin(3) - t64;
t48 = qJD(1) * t120 + qJD(5) + t59;
t23 = -t67 * pkin(5) + t69 * pkin(8) + t48;
t4 = t129 * t11 + t215 * t12;
t182 = qJDD(4) * pkin(8) + qJD(6) * t23 + t4;
t201 = qJD(2) * t131;
t169 = -qJD(5) + t201;
t79 = t131 * qJ(2) + t130 * t138;
t76 = -pkin(7) + t79;
t216 = qJ(5) - t76;
t176 = qJD(4) * t216;
t143 = -t169 * t135 + t137 * t176;
t40 = t135 * t176 + t169 * t137;
t19 = t129 * t143 + t215 * t40;
t229 = sin(qJ(1));
t230 = cos(qJ(1));
t72 = -t229 * t130 - t230 * t131;
t232 = g(1) * t72;
t61 = t216 * t137;
t25 = t216 * t211 - t215 * t61;
t78 = -t130 * qJ(2) + t131 * t138;
t75 = pkin(3) - t78;
t154 = t120 + t75;
t27 = pkin(5) * t152 + t73 * pkin(8) + t154;
t234 = -t1 * t73 - t13 * t70 - (qJD(6) * t27 + t19) * t210 - t182 * t152 + t25 * t43 - t232;
t100 = t129 * pkin(4) + pkin(8);
t125 = qJ(4) + pkin(10);
t109 = sin(t125);
t110 = cos(t125);
t74 = t230 * t130 - t229 * t131;
t167 = g(2) * t74 + t232;
t233 = t109 * t167 + (-pkin(4) * t202 - t69 * pkin(5) - t67 * pkin(8) + qJD(6) * t100) * t210 - g(3) * t110 + t1;
t231 = g(2) * t72;
t228 = g(3) * t137;
t53 = t134 * qJD(4) - t136 * t69;
t227 = t53 * t69;
t226 = t69 * t51;
t32 = t215 * t42;
t16 = t129 * t35 + t32;
t225 = t130 * t70 - t131 * t69;
t224 = -t130 * t68 - t152 * t203;
t222 = t110 * t72;
t221 = t110 * t74;
t199 = qJD(6) * t134;
t195 = qJD(1) * qJD(4);
t183 = t135 * t195;
t45 = qJD(4) * t89 + t73 * qJDD(1) - t129 * t183;
t21 = qJD(6) * t198 + t134 * qJDD(4) - t136 * t45 + t69 * t199;
t217 = t21 * t134;
t214 = pkin(1) * qJDD(1);
t14 = qJD(4) * pkin(8) + t16;
t213 = qJD(6) * t14;
t209 = qJDD(3) + g(3);
t208 = t230 * pkin(1) + t229 * qJ(2);
t207 = g(1) * t229 - g(2) * t230;
t126 = t135 ^ 2;
t206 = -t137 ^ 2 + t126;
t139 = qJD(4) ^ 2;
t140 = qJD(1) ^ 2;
t205 = t139 + t140;
t197 = qJ(2) * qJDD(1);
t194 = qJDD(4) * t135;
t193 = qJDD(4) * t137;
t192 = t130 * qJDD(1);
t189 = t137 * qJDD(1);
t188 = 0.2e1 * t196;
t187 = t73 * t199;
t186 = t230 * pkin(2) + t208;
t91 = t130 * t196;
t181 = -qJ(2) * t192 + t131 * t85;
t177 = -t136 * qJDD(4) - t134 * t45;
t171 = 0.2e1 * t137 * t195;
t170 = qJDD(2) - t214;
t49 = t181 - t91;
t168 = g(1) * t74 - t231;
t166 = -t229 * pkin(1) + t230 * qJ(2);
t165 = -pkin(4) * t200 + t130 * qJD(2);
t164 = -qJD(6) * t131 + t224;
t162 = t152 * t21 - t68 * t53;
t22 = t53 * qJD(6) + t177;
t161 = -t152 * t22 + t68 * t51;
t159 = -t213 + t231;
t158 = t64 * t130 - t65 * t131;
t46 = qJDD(1) * pkin(3) - t49;
t63 = t152 * t130;
t156 = qJD(6) * t63 + t204;
t155 = -t136 * t43 + (t134 * t67 - t199) * t210;
t153 = -g(1) * t221 - t27 * t43;
t151 = g(1) * t230 + g(2) * t229;
t150 = t59 * qJD(1) - t167 - t47;
t149 = -t229 * pkin(2) + t166;
t20 = t215 * t41 - t220;
t148 = t100 * t43 + (t13 + t20) * t210;
t147 = -g(2) * t221 - g(3) * t109 - t182;
t144 = -qJDD(4) * t76 + (-qJD(1) * t75 - t201 - t59) * qJD(4);
t142 = t160 * t136 + t187 * t210;
t28 = qJDD(5) + t46 + (-t183 + t189) * pkin(4);
t141 = qJDD(1) * t75 - t139 * t76 - t168 + t46 + t91;
t133 = -qJ(5) - pkin(7);
t105 = t120 + pkin(3);
t101 = -t215 * pkin(4) - pkin(5);
t83 = -t139 * t135 + t193;
t82 = -t139 * t137 - t194;
t62 = t73 * t130;
t37 = t74 * t134 - t136 * t222;
t36 = t134 * t222 + t74 * t136;
t26 = -t68 * pkin(5) + t70 * pkin(8) + t165;
t24 = -t129 * t61 - t216 * t180;
t18 = t129 * t41 + t32;
t17 = t129 * t40 - t215 * t143;
t8 = -t44 * pkin(5) + t45 * pkin(8) + t28;
t7 = t136 * t8;
t6 = t134 * t23 + t136 * t14;
t5 = -t134 * t14 + t136 * t23;
t2 = [qJDD(1), t207, t151, -qJDD(2) + t207 + 0.2e1 * t214, -t151 + t188 + 0.2e1 * t197, -t170 * pkin(1) - g(1) * t166 - g(2) * t208 + (t188 + t197) * qJ(2), -t78 * qJDD(1) - t168 - t181 + 0.2e1 * t91, t79 * qJDD(1) + t167 + t223 + 0.2e1 * t93, -g(1) * t149 - g(2) * t186 - t158 * qJD(2) + t49 * t78 + t50 * t79, t126 * qJDD(1) + t135 * t171, 0.2e1 * t135 * t189 - 0.2e1 * t206 * t195, t82, -t83, 0, t144 * t135 + t141 * t137, -t141 * t135 + t144 * t137, t15 * t70 - t152 * t4 + t16 * t68 - t17 * t69 + t19 * t67 - t24 * t45 + t25 * t44 + t3 * t73 - t167, t4 * t25 + t16 * t19 - t3 * t24 - t15 * t17 + t28 * t154 + t48 * t165 - g(1) * (t74 * t105 - t72 * t133 + t149) - g(2) * (-t72 * t105 - t74 * t133 + t186) t53 * t187 + (-t21 * t73 - t53 * t70) * t136 (t134 * t53 + t136 * t51) * t70 + (t217 + t136 * t22 + (-t134 * t51 + t136 * t53) * qJD(6)) * t73, t142 + t162, t161 + t235, -t152 * t43 - t210 * t68, -g(2) * t37 + t17 * t51 + t24 * t22 - t5 * t68 + t7 * t152 + (t26 * t210 + (-t13 * t73 - t14 * t152 - t210 * t25) * qJD(6) + t153) * t136 + t234 * t134, -g(2) * t36 + t17 * t53 + t24 * t21 + t6 * t68 + (-(-qJD(6) * t25 + t26) * t210 - (t8 - t213) * t152 + t13 * t212 - t153) * t134 + t234 * t136; 0, 0, 0, -qJDD(1), -t140, -t140 * qJ(2) + t170 - t207, -t130 * t140 - t191, -t131 * t140 + t192, t158 * qJD(1) + t50 * t130 + t49 * t131 - t207, 0, 0, 0, 0, 0 (0.2e1 * t183 - t189) * t131 + (-t205 * t137 - t194) * t130 (t171 + t190) * t131 + (t205 * t135 - t193) * t130, t224 * t67 - t225 * t69 + t63 * t44 - t62 * t45, -t28 * t131 - t225 * t15 + t224 * t16 - t48 * t204 - t3 * t62 + t4 * t63 - t207, 0, 0, 0, 0, 0 -(-t136 * t131 - t134 * t63) * t43 + t62 * t22 - (t134 * t164 + t136 * t156) * t210 + t225 * t51 (-t134 * t131 + t136 * t63) * t43 + t62 * t21 - (-t134 * t156 + t136 * t164) * t210 + t225 * t53; 0, 0, 0, 0, 0, 0, 0, 0, t209, 0, 0, 0, 0, 0, t83, t82, t152 * t45 + t73 * t44 + t70 * t67 - t68 * t69, -t15 * t68 + t152 * t3 + t16 * t70 + t4 * t73 + g(3), 0, 0, 0, 0, 0, t161 - t235, t142 - t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135 * t140 * t137, t206 * t140, -t190, -t189, qJDD(4), t150 * t135 + t115 + t228, -t209 * t135 + t150 * t137 (-t16 + t18) * t69 + (t15 - t20) * t67 + (t129 * t44 + t215 * t45) * pkin(4), t15 * t18 - t16 * t20 + (t215 * t3 + t228 + t129 * t4 + (qJD(1) * t48 - t167) * t135) * pkin(4), t236 * t53 + t217 (t21 - t237) * t136 + (-t210 * t53 - t22) * t134, t227 - t238, t155 - t226, t210 * t69, t101 * t22 + t148 * t134 - t136 * t233 - t18 * t51 + t5 * t69, t101 * t21 + t134 * t233 + t148 * t136 - t18 * t53 - t6 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 ^ 2 - t69 ^ 2, -t15 * t69 - t16 * t67 - t168 + t28, 0, 0, 0, 0, 0, t155 + t226, t227 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t51, -t51 ^ 2 + t53 ^ 2, t21 + t237, -t177 + (-qJD(6) + t210) * t53, -t43, -g(1) * t36 - t13 * t53 + t134 * t147 + t136 * t159 + t210 * t6 + t7, g(1) * t37 + t13 * t51 + t5 * t210 + (-t159 - t8) * t134 + t147 * t136;];
tau_reg  = t2;
