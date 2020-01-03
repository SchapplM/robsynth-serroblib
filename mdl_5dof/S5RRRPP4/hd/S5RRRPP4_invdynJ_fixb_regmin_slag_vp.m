% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:48
% EndTime: 2019-12-31 20:55:54
% DurationCPUTime: 1.98s
% Computational Cost: add. (3705->302), mult. (8795->382), div. (0->0), fcn. (6167->12), ass. (0->172)
t153 = qJ(2) + qJ(3);
t144 = sin(t153);
t145 = cos(t153);
t158 = sin(qJ(1));
t161 = cos(qJ(1));
t190 = g(1) * t161 + g(2) * t158;
t248 = -g(3) * t145 + t144 * t190;
t150 = qJD(2) + qJD(3);
t154 = sin(pkin(8));
t155 = cos(pkin(8));
t159 = cos(qJ(3));
t160 = cos(qJ(2));
t216 = qJD(1) * t160;
t203 = t159 * t216;
t156 = sin(qJ(3));
t157 = sin(qJ(2));
t217 = qJD(1) * t157;
t204 = t156 * t217;
t81 = -t203 + t204;
t83 = -t156 * t216 - t159 * t217;
t55 = t154 * t83 - t155 * t81;
t247 = t55 * t150;
t213 = qJD(1) * qJD(2);
t200 = t160 * t213;
t212 = t157 * qJDD(1);
t245 = t200 + t212;
t183 = -t154 * t81 - t155 * t83;
t244 = t183 ^ 2;
t241 = pkin(6) + pkin(7);
t146 = t160 * pkin(2);
t231 = pkin(1) + t146;
t111 = t241 * t160;
t101 = qJD(1) * t111;
t88 = t159 * t101;
t110 = t241 * t157;
t99 = qJD(1) * t110;
t196 = t156 * t99 - t88;
t223 = t81 * qJ(4);
t181 = t196 + t223;
t210 = t154 * t156 * pkin(2);
t214 = qJD(3) * t159;
t84 = t156 * t101;
t227 = -t159 * t99 - t84;
t77 = t83 * qJ(4);
t49 = t77 + t227;
t228 = -qJD(3) * t210 - t154 * t181 + (pkin(2) * t214 - t49) * t155;
t225 = qJD(2) * pkin(2);
t91 = -t99 + t225;
t194 = t159 * t91 - t84;
t45 = t194 + t77;
t143 = pkin(8) + t153;
t130 = sin(t143);
t131 = cos(t143);
t185 = t131 * pkin(4) + t130 * qJ(5);
t243 = g(1) * t158 - g(2) * t161;
t148 = qJDD(2) + qJDD(3);
t182 = -t156 * t91 - t88;
t64 = qJDD(2) * pkin(2) - t241 * t245;
t201 = t157 * t213;
t211 = t160 * qJDD(1);
t66 = t241 * (-t201 + t211);
t171 = qJD(3) * t182 - t156 * t66 + t159 * t64;
t47 = qJD(3) * t203 - t150 * t204 + t156 * t211 + t245 * t159;
t13 = t148 * pkin(3) - t47 * qJ(4) + t83 * qJD(4) + t171;
t215 = qJD(3) * t156;
t242 = (qJD(3) * t91 + t66) * t159 - t101 * t215 + t156 * t64;
t184 = t156 * t212 - t159 * t211;
t94 = t156 * t160 + t159 * t157;
t63 = t150 * t94;
t48 = qJD(1) * t63 + t184;
t16 = -t48 * qJ(4) - t81 * qJD(4) + t242;
t3 = t155 * t13 - t154 * t16;
t2 = -t148 * pkin(4) + qJDD(5) - t3;
t109 = t231 * qJD(1);
t65 = t81 * pkin(3) + qJD(4) - t109;
t26 = -pkin(4) * t55 - qJ(5) * t183 + t65;
t169 = -g(3) * t131 + t190 * t130 - t183 * t26 - t2;
t240 = t83 * pkin(3);
t239 = pkin(3) * t144;
t238 = pkin(4) * t130;
t46 = -t182 - t223;
t41 = t155 * t46;
t21 = t154 * t45 + t41;
t233 = t21 * t183;
t232 = t83 * t81;
t4 = t154 * t13 + t155 * t16;
t40 = t150 * pkin(3) + t45;
t20 = t154 * t40 + t41;
t230 = qJD(5) + t228;
t221 = t155 * t156;
t229 = -t154 * t49 + t155 * t181 + (t154 * t159 + t221) * qJD(3) * pkin(2);
t226 = -t156 * t110 + t159 * t111;
t224 = t154 * t46;
t222 = qJ(5) * t131;
t22 = t155 * t45 - t224;
t220 = qJD(5) - t22;
t136 = t159 * pkin(2) + pkin(3);
t76 = pkin(2) * t221 + t154 * t136;
t134 = pkin(3) * t145;
t219 = t134 + t146;
t151 = t157 ^ 2;
t218 = -t160 ^ 2 + t151;
t209 = t148 * qJ(5) + t4;
t140 = t157 * t225;
t208 = t134 + t185;
t206 = qJD(2) * t241;
t205 = t229 * t183;
t202 = t63 * pkin(3) + t140;
t103 = -t157 * pkin(2) - t239;
t199 = t103 - t238;
t23 = t154 * t47 + t155 * t48;
t193 = -t159 * t110 - t156 * t111;
t93 = t156 * t157 - t159 * t160;
t192 = t93 * pkin(3) - t231;
t191 = -t238 - t239;
t19 = t155 * t40 - t224;
t17 = -t150 * pkin(4) + qJD(5) - t19;
t18 = t150 * qJ(5) + t20;
t188 = -t17 * t55 + t18 * t183;
t187 = t183 * t20 + t19 * t55;
t186 = -t55 ^ 2 - t244;
t24 = -t154 * t48 + t155 * t47;
t75 = t155 * t136 - t210;
t180 = -t94 * qJ(4) + t193;
t179 = -0.2e1 * pkin(1) * t213 - pkin(6) * qJDD(2);
t100 = t157 * t206;
t102 = t160 * t206;
t178 = -t159 * t100 - t156 * t102 - t110 * t214 - t111 * t215;
t78 = pkin(2) * t201 - qJDD(1) * t231;
t31 = pkin(4) * t183 - t55 * qJ(5) - t240;
t162 = qJD(2) ^ 2;
t174 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t162 + t243;
t163 = qJD(1) ^ 2;
t173 = pkin(1) * t163 - pkin(6) * qJDD(1) + t190;
t172 = t48 * pkin(3) + qJDD(4) + t78;
t52 = -t93 * qJ(4) + t226;
t33 = t154 * t52 - t155 * t180;
t34 = t154 * t180 + t155 * t52;
t168 = -qJD(3) * t226 + t156 * t100 - t159 * t102;
t62 = t150 * t93;
t165 = t62 * qJ(4) - t94 * qJD(4) + t168;
t27 = -t63 * qJ(4) - t93 * qJD(4) + t178;
t7 = t154 * t27 - t155 * t165;
t8 = t154 * t165 + t155 * t27;
t170 = t183 * t7 - t34 * t23 + t33 * t24 + t55 * t8 - t190;
t167 = -g(3) * t130 - t131 * t190 + t26 * t55 + t209;
t166 = g(3) * t144 - t109 * t81 + t190 * t145 - t242;
t5 = t23 * pkin(4) - t24 * qJ(5) - qJD(5) * t183 + t172;
t164 = -t109 * t83 + t171 + t248;
t149 = -qJ(4) - t241;
t141 = t150 * qJD(5);
t139 = pkin(2) * t217;
t132 = -t155 * pkin(3) - pkin(4);
t129 = t154 * pkin(3) + qJ(5);
t105 = t161 * t222;
t104 = t158 * t222;
t98 = pkin(1) + t219;
t90 = t161 * t98;
t71 = -pkin(4) - t75;
t70 = qJ(5) + t76;
t60 = -t154 * t93 + t155 * t94;
t59 = t154 * t94 + t155 * t93;
t50 = -t81 ^ 2 + t83 ^ 2;
t38 = -t154 * t63 - t155 * t62;
t37 = -t154 * t62 + t155 * t63;
t36 = -t184 + (-qJD(1) * t94 - t83) * t150;
t35 = t81 * t150 + t47;
t32 = t59 * pkin(4) - t60 * qJ(5) + t192;
t30 = t139 + t31;
t9 = t37 * pkin(4) - t38 * qJ(5) - t60 * qJD(5) + t202;
t1 = t141 + t209;
t6 = [qJDD(1), t243, t190, t151 * qJDD(1) + 0.2e1 * t157 * t200, 0.2e1 * t157 * t211 - 0.2e1 * t213 * t218, qJDD(2) * t157 + t162 * t160, qJDD(2) * t160 - t162 * t157, 0, t157 * t179 + t160 * t174, -t157 * t174 + t160 * t179, t47 * t94 + t83 * t62, -t47 * t93 - t94 * t48 + t62 * t81 + t83 * t63, t94 * t148 - t62 * t150, -t93 * t148 - t63 * t150, 0, -t109 * t63 + t140 * t81 + t145 * t243 + t148 * t193 + t150 * t168 - t231 * t48 + t78 * t93, t109 * t62 - t140 * t83 - t144 * t243 - t148 * t226 - t150 * t178 - t231 * t47 + t78 * t94, -t19 * t38 - t20 * t37 - t3 * t60 - t4 * t59 + t170, t4 * t34 + t20 * t8 - t3 * t33 - t19 * t7 + t172 * t192 + t65 * t202 - g(1) * (-t161 * t149 - t158 * t98) - g(2) * (-t158 * t149 + t90), t131 * t243 - t33 * t148 - t7 * t150 + t32 * t23 + t26 * t37 + t5 * t59 - t55 * t9, -t1 * t59 + t17 * t38 - t18 * t37 + t2 * t60 + t170, t130 * t243 + t34 * t148 + t8 * t150 - t183 * t9 - t32 * t24 - t26 * t38 - t5 * t60, -g(2) * t90 + t1 * t34 + t17 * t7 + t18 * t8 + t2 * t33 + t26 * t9 + t5 * t32 + (g(1) * t149 - g(2) * t185) * t161 + (-g(1) * (-t185 - t98) + g(2) * t149) * t158; 0, 0, 0, -t157 * t163 * t160, t218 * t163, t212, t211, qJDD(2), -g(3) * t160 + t157 * t173, g(3) * t157 + t160 * t173, -t232, t50, t35, t36, t148, -t196 * t150 + (t159 * t148 - t150 * t215 - t217 * t81) * pkin(2) + t164, t227 * t150 + (-t156 * t148 - t150 * t214 + t217 * t83) * pkin(2) + t166, t228 * t55 - t76 * t23 - t75 * t24 + t187 + t205, t4 * t76 + t3 * t75 - t65 * (t139 - t240) - g(3) * t219 + t228 * t20 - t229 * t19 - t190 * t103, -t71 * t148 - t229 * t150 + t30 * t55 + t169, -t70 * t23 + t230 * t55 + t71 * t24 + t188 + t205, t70 * t148 + t230 * t150 + t183 * t30 + t141 + t167, t1 * t70 + t2 * t71 - t26 * t30 - g(1) * (t161 * t199 + t105) - g(2) * (t158 * t199 + t104) - g(3) * (t146 + t208) + t230 * t18 + t229 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t232, t50, t35, t36, t148, -t150 * t182 + t164, t150 * t194 + t166, -t233 - t22 * t55 + (-t154 * t23 - t155 * t24) * pkin(3) + t187, t19 * t21 - t20 * t22 + (t154 * t4 + t155 * t3 + t65 * t83 + t248) * pkin(3), -t132 * t148 + t21 * t150 + t31 * t55 + t169, -t129 * t23 + t132 * t24 + t220 * t55 + t188 - t233, t129 * t148 - t22 * t150 + t183 * t31 + 0.2e1 * t141 + t167, t1 * t129 + t2 * t132 - t26 * t31 - t17 * t21 - g(1) * (t161 * t191 + t105) - g(2) * (t158 * t191 + t104) - g(3) * t208 + t220 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, t183 * t19 - t20 * t55 + t172 - t243, t150 * t183 + t23, t186, -t24 - t247, -t17 * t183 - t18 * t55 - t243 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183 * t55 - t148, t24 - t247, -t150 ^ 2 - t244, -t18 * t150 - t169;];
tau_reg = t6;
