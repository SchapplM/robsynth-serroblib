% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:24
% EndTime: 2019-03-09 01:49:30
% DurationCPUTime: 2.27s
% Computational Cost: add. (2106->372), mult. (3875->477), div. (0->0), fcn. (2488->10), ass. (0->189)
t130 = sin(qJ(4));
t208 = qJD(1) * t130;
t100 = qJD(6) + t208;
t129 = sin(qJ(6));
t132 = cos(qJ(6));
t125 = sin(pkin(9));
t133 = cos(qJ(4));
t207 = qJD(1) * t133;
t182 = t125 * t207;
t126 = cos(pkin(9));
t200 = t126 * qJD(4);
t73 = t182 - t200;
t181 = t126 * t207;
t205 = qJD(4) * t125;
t75 = t181 + t205;
t152 = t129 * t73 - t132 * t75;
t246 = t100 * t152;
t131 = sin(qJ(1));
t134 = cos(qJ(1));
t242 = g(1) * t131 - g(2) * t134;
t222 = qJ(5) * t133;
t157 = pkin(4) * t130 - t222;
t118 = qJDD(1) * qJ(3);
t124 = qJDD(1) * pkin(1);
t189 = t124 - qJDD(2);
t173 = -t118 - t189;
t158 = pkin(4) * t133 + qJ(5) * t130;
t60 = t158 * qJD(4) - t133 * qJD(5) + qJD(3);
t18 = t60 * qJD(1) + t157 * qJDD(1) - t173;
t103 = qJ(2) * qJD(1) + qJD(3);
t94 = -pkin(7) * qJD(1) + t103;
t226 = t133 * t94;
t119 = qJD(1) * qJD(2);
t120 = qJ(2) * qJDD(1);
t172 = qJDD(3) + t119 + t120;
t81 = -pkin(7) * qJDD(1) + t172;
t27 = qJDD(4) * qJ(5) + t130 * t81 + (qJD(5) + t226) * qJD(4);
t8 = -t125 * t27 + t126 * t18;
t9 = t125 * t18 + t126 * t27;
t245 = -t125 * t9 - t126 * t8 - t242;
t127 = -pkin(7) + qJ(2);
t203 = qJD(4) * t133;
t206 = qJD(2) * t130;
t244 = t127 * t203 + t206;
t214 = t132 * t126;
t78 = t125 * t129 - t214;
t145 = t100 * t78;
t165 = g(1) * t134 + g(2) * t131;
t243 = t165 * t126;
t79 = t125 * t132 + t126 * t129;
t143 = t79 * qJD(6);
t201 = qJD(6) * t132;
t202 = qJD(6) * t129;
t241 = -t125 * t202 + t126 * t201;
t128 = pkin(1) + qJ(3);
t240 = qJD(1) * t128;
t239 = -qJD(6) + t100;
t142 = t79 * qJD(1);
t229 = t130 * t142 + t143;
t198 = qJD(1) * qJD(4);
t175 = t133 * t198;
t191 = t130 * qJDD(1);
t144 = t175 + t191;
t77 = qJDD(6) + t144;
t233 = t78 * t77;
t238 = -t229 * t100 - t233;
t235 = g(3) * t130;
t237 = t165 * t133 - t235;
t95 = -qJD(2) + t240;
t236 = qJD(4) * (qJD(2) + t95 + t240) + qJDD(4) * t127;
t190 = t133 * qJDD(1);
t168 = -t126 * qJDD(4) + t125 * t190;
t176 = t130 * t198;
t44 = t125 * t176 - t168;
t223 = t125 * qJDD(4) + t126 * t190;
t45 = t126 * t176 - t223;
t7 = -qJD(6) * t152 - t129 * t45 - t132 * t44;
t112 = 0.2e1 * t119;
t234 = g(3) * t133;
t232 = t79 * t77;
t231 = pkin(8) + qJ(5);
t85 = t157 + t128;
t52 = qJD(1) * t85 - qJD(2);
t84 = t130 * t94;
t68 = qJD(4) * qJ(5) + t84;
t20 = t125 * t52 + t126 * t68;
t183 = t125 * t208;
t230 = -t129 * t183 + t208 * t214 + t241;
t219 = t126 * t133;
t83 = t158 * qJD(1);
t36 = t125 * t83 + t94 * t219;
t217 = t127 * t130;
t41 = t125 * t85 + t126 * t217;
t30 = t129 * t75 + t132 * t73;
t228 = t100 * t30;
t204 = qJD(4) * t130;
t160 = -qJDD(4) * pkin(4) + t94 * t204 + qJDD(5);
t29 = -t133 * t81 + t160;
t227 = t133 * t29;
t225 = t29 * t125;
t224 = t29 * t126;
t221 = t125 * t133;
t136 = qJD(1) ^ 2;
t220 = t125 * t136;
t218 = t126 * t136;
t216 = t130 * t131;
t215 = t130 * t134;
t62 = -qJD(4) * pkin(4) + qJD(5) - t226;
t213 = -qJD(5) + t62;
t212 = t134 * pkin(1) + t131 * qJ(2);
t122 = t130 ^ 2;
t123 = t133 ^ 2;
t210 = t122 - t123;
t135 = qJD(4) ^ 2;
t209 = -t135 - t136;
t199 = t133 * qJD(2);
t197 = qJD(3) * qJD(1);
t196 = qJDD(1) * t125;
t195 = qJDD(1) * t126;
t194 = qJDD(1) * t128;
t192 = qJDD(4) * t130;
t24 = t125 * t60 + t126 * t244;
t188 = pkin(8) * t126 * t130;
t187 = t30 * t207;
t186 = t152 * t207;
t185 = t134 * qJ(3) + t212;
t4 = pkin(5) * t144 + pkin(8) * t45 + t8;
t5 = pkin(8) * t44 + t9;
t184 = -t129 * t5 + t132 * t4;
t180 = t125 * t204;
t174 = qJDD(2) - t242;
t171 = -t125 * t127 + pkin(5);
t170 = pkin(5) * t125 - t127;
t19 = -t125 * t68 + t126 * t52;
t167 = -0.2e1 * t175;
t166 = -t124 + t174;
t35 = t126 * t83 - t94 * t221;
t162 = -t8 * t125 + t9 * t126;
t161 = t129 * t4 + t132 * t5;
t10 = pkin(5) * t208 - pkin(8) * t75 + t19;
t12 = -pkin(8) * t73 + t20;
t1 = t10 * t132 - t12 * t129;
t2 = t10 * t129 + t12 * t132;
t156 = t125 * t20 + t126 * t19;
t155 = t125 * t19 - t126 * t20;
t70 = t126 * t85;
t28 = -pkin(8) * t219 + t171 * t130 + t70;
t33 = -pkin(8) * t221 + t41;
t154 = -t129 * t33 + t132 * t28;
t153 = t129 * t28 + t132 * t33;
t151 = -t118 + t166;
t150 = t165 * t125;
t149 = t230 * t100 + t232;
t92 = t231 * t126;
t148 = qJD(5) * t125 + qJD(6) * t92 + (pkin(5) * t133 + t188) * qJD(1) + t35;
t91 = t231 * t125;
t147 = pkin(8) * t183 - qJD(5) * t126 + qJD(6) * t91 + t36;
t6 = t129 * t44 - t132 * t45 - t73 * t201 - t75 * t202;
t146 = t112 + 0.2e1 * t120 - t165;
t141 = qJD(1) * t95 + t165;
t140 = t141 - t81;
t138 = -t165 * t130 - t234;
t82 = -t173 + t197;
t137 = -t127 * t135 + t194 + t197 + t242 + t82;
t117 = pkin(9) + qJ(6);
t111 = t134 * qJ(2);
t108 = qJDD(4) * t133;
t107 = cos(t117);
t106 = sin(t117);
t102 = -pkin(5) * t126 - pkin(4);
t72 = t170 * t133;
t59 = t78 * t133;
t58 = t79 * t133;
t56 = -t106 * t131 + t107 * t215;
t55 = -t106 * t215 - t107 * t131;
t54 = -t106 * t134 - t107 * t216;
t53 = t106 * t216 - t107 * t134;
t51 = -pkin(5) * t183 + t84;
t48 = -t170 * t204 - t199;
t47 = t126 * t60;
t40 = -t125 * t217 + t70;
t34 = pkin(5) * t73 + t62;
t23 = -t125 * t244 + t47;
t22 = -t129 * t130 * t200 - t132 * t180 + t241 * t133;
t21 = -t133 * t143 + t78 * t204;
t14 = pkin(8) * t180 + t24;
t13 = -t125 * t206 + t47 + (t171 * t133 + t188) * qJD(4);
t11 = -pkin(5) * t44 + t29;
t3 = [qJDD(1), t242, t165, -0.2e1 * t124 + t174, t146, t189 * pkin(1) - g(1) * (-pkin(1) * t131 + t111) - g(2) * t212 + (t112 + t120) * qJ(2), qJDD(3) + t146, -t151 + t194 + 0.2e1 * t197, t82 * t128 + t95 * qJD(3) + t172 * qJ(2) + t103 * qJD(2) - g(1) * (-t128 * t131 + t111) - g(2) * t185, qJDD(1) * t123 + t130 * t167, -0.2e1 * t130 * t190 + 0.2e1 * t210 * t198, -t130 * t135 + t108, -t133 * t135 - t192, 0, t137 * t130 + t236 * t133, -t236 * t130 + t137 * t133, t150 + (-qJD(2) * t73 + t225 + t127 * t44 + (qJD(1) * t40 + t19) * qJD(4)) * t133 + (t23 * qJD(1) + t40 * qJDD(1) + t8 + t242 * t126 + (-t125 * t62 + t127 * t73) * qJD(4)) * t130, t243 + (-qJD(2) * t75 + t224 + t127 * t45 + (-qJD(1) * t41 - t20) * qJD(4)) * t133 + (-t24 * qJD(1) - t41 * qJDD(1) - t9 - t242 * t125 + (-t126 * t62 + t127 * t75) * qJD(4)) * t130, t133 * t245 + t156 * t204 - t23 * t75 - t24 * t73 + t40 * t45 + t41 * t44, t9 * t41 + t20 * t24 + t8 * t40 + t19 * t23 - t62 * t199 - g(1) * (-pkin(7) * t134 + t111) - g(2) * (pkin(4) * t215 - t134 * t222 + t185) + (t62 * t204 - t227) * t127 + (g(2) * pkin(7) + g(1) * t85) * t131, -t152 * t21 - t59 * t6, t152 * t22 - t21 * t30 - t58 * t6 + t59 * t7, t100 * t21 + t130 * t6 - t152 * t203 - t59 * t77, -t100 * t22 - t130 * t7 - t30 * t203 - t58 * t77, t100 * t203 + t130 * t77 (-t129 * t14 + t132 * t13) * t100 + t154 * t77 + t184 * t130 + t1 * t203 + t48 * t30 + t72 * t7 + t11 * t58 + t34 * t22 - g(1) * t54 - g(2) * t56 + (-t100 * t153 - t130 * t2) * qJD(6) -(t129 * t13 + t132 * t14) * t100 - t153 * t77 - t161 * t130 - t2 * t203 - t48 * t152 + t72 * t6 - t11 * t59 + t34 * t21 - g(1) * t53 - g(2) * t55 + (-t1 * t130 - t100 * t154) * qJD(6); 0, 0, 0, qJDD(1), -t136, -qJ(2) * t136 + t166, -t136, -qJDD(1) (-qJD(3) - t103) * qJD(1) + t151, 0, 0, 0, 0, 0, t167 - t191, 0.2e1 * t176 - t190, -t126 * t191 + t122 * t220 + (t73 - t200) * t207, t125 * t191 + t122 * t218 + (t75 + t205) * t207, -t125 * t44 - t126 * t45 + (-t125 * t75 + t126 * t73) * t208 (t130 * t155 + t133 * t62) * qJD(1) + t245, 0, 0, 0, 0, 0, t187 - t238, t149 - t186; 0, 0, 0, 0, 0, 0, qJDD(1), -t136, -t141 + t172, 0, 0, 0, 0, 0, t209 * t130 + t108, t209 * t133 - t192, -t122 * t196 + t133 * t44 + (-t218 + (t73 - 0.2e1 * t182) * qJD(4)) * t130, -t122 * t195 + t133 * t45 + (t220 + (t75 - 0.2e1 * t181) * qJD(4)) * t130 (qJD(1) * t75 + t130 * t44 - t73 * t203) * t126 + (qJD(1) * t73 - t130 * t45 + t75 * t203) * t125, -t227 + t162 * t130 - t156 * qJD(1) + (t130 * t62 - t133 * t155) * qJD(4) - t165, 0, 0, 0, 0, 0, qJD(1) * t145 + (-qJD(4) * t100 * t79 - t7) * t133 + (qJD(4) * t30 + qJD(6) * t145 - t232) * t130, t100 * t142 + (qJD(4) * t145 - t6) * t133 + (-qJD(4) * t152 + t100 * t143 + t233) * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 * t136 * t130, -t210 * t136, t190, -t191, qJDD(4), -t133 * t140 + t235, t130 * t140 + t234, pkin(4) * t44 - t224 + (-t243 + (-qJ(5) * t205 - t19) * qJD(1)) * t133 + (-qJ(5) * t196 + g(3) * t126 - t73 * t94 + (t213 * t125 - t35) * qJD(1)) * t130, pkin(4) * t45 + t225 + (t150 + (-qJ(5) * t200 + t20) * qJD(1)) * t133 + (-qJ(5) * t195 - g(3) * t125 - t75 * t94 + (t213 * t126 + t36) * qJD(1)) * t130, t35 * t75 + t36 * t73 + (qJ(5) * t44 - qJD(5) * t73 - t19 * t208 + t9) * t126 + (-qJ(5) * t45 + qJD(5) * t75 - t20 * t208 - t8) * t125 + t138, -t62 * t84 - t19 * t35 - t20 * t36 - t155 * qJD(5) + (-t29 - t237) * pkin(4) + (t138 + t162) * qJ(5), -t152 * t230 + t6 * t79, t152 * t229 - t230 * t30 - t6 * t78 - t7 * t79, t149 + t186, t187 + t238, -t100 * t207 (-t129 * t92 - t132 * t91) * t77 + t102 * t7 + t11 * t78 - t1 * t207 - t51 * t30 + t229 * t34 + (t129 * t147 - t132 * t148) * t100 - t237 * t107 -(-t129 * t91 + t132 * t92) * t77 + t102 * t6 + t11 * t79 + t2 * t207 + t51 * t152 + t230 * t34 + (t129 * t148 + t132 * t147) * t100 + t237 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t75 - t205) * t208 + t168 (-t73 - t200) * t208 + t223, -t73 ^ 2 - t75 ^ 2, -t235 + t19 * t75 + t20 * t73 + (t165 - t81) * t133 + t160, 0, 0, 0, 0, 0, t7 - t246, t6 - t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 * t30, t152 ^ 2 - t30 ^ 2, t6 + t228, -t7 - t246, t77, -g(1) * t55 + g(2) * t53 + t106 * t234 + t152 * t34 + t239 * t2 + t184, g(1) * t56 - g(2) * t54 + t239 * t1 + t107 * t234 + t30 * t34 - t161;];
tau_reg  = t3;
