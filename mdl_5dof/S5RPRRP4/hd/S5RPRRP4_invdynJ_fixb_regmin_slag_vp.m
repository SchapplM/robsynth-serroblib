% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:27
% EndTime: 2020-01-03 11:50:34
% DurationCPUTime: 2.29s
% Computational Cost: add. (2117->290), mult. (5124->395), div. (0->0), fcn. (3628->10), ass. (0->182)
t129 = sin(pkin(8));
t131 = sin(qJ(4));
t183 = t129 * qJDD(1);
t132 = sin(qJ(3));
t185 = qJDD(1) * t132;
t134 = cos(qJ(4));
t135 = cos(qJ(3));
t205 = t134 * t135;
t181 = qJD(3) + qJD(4);
t79 = t131 * t135 + t134 * t132;
t244 = t181 * t79;
t19 = (qJD(1) * t244 + t131 * t185) * t129 - t183 * t205;
t124 = t129 ^ 2;
t237 = 0.2e1 * t124;
t216 = qJDD(1) * pkin(1);
t112 = qJDD(2) - t216;
t133 = sin(qJ(1));
t136 = cos(qJ(1));
t240 = -g(2) * t136 - g(3) * t133;
t243 = t240 - t112;
t130 = cos(pkin(8));
t190 = t130 * qJD(1);
t105 = -qJD(3) + t190;
t212 = t129 * t135;
t180 = pkin(7) * t212;
t219 = qJ(2) * t132;
t148 = -t130 * t219 - t180;
t85 = -t130 * pkin(2) - t129 * pkin(6) - pkin(1);
t67 = t85 * qJD(1) + qJD(2);
t59 = t135 * t67;
t38 = t148 * qJD(1) + t59;
t24 = -t105 * pkin(3) + t38;
t218 = qJ(2) * t135;
t104 = t130 * t218;
t196 = qJD(1) * t129;
t177 = t132 * t196;
t39 = -pkin(7) * t177 + qJD(1) * t104 + t132 * t67;
t31 = t131 * t39;
t169 = t134 * t24 - t31;
t158 = t131 * t177;
t195 = qJD(1) * t135;
t176 = t129 * t195;
t54 = t134 * t176 - t158;
t48 = t54 * qJ(5);
t242 = t48 - t169;
t186 = qJDD(1) * qJ(2);
t188 = qJD(1) * qJD(2);
t151 = t186 + t188;
t241 = t130 * t151;
t128 = qJ(3) + qJ(4);
t115 = sin(t128);
t232 = g(1) * t129;
t116 = cos(t128);
t203 = t136 * t116;
t209 = t133 * t115;
t60 = -t130 * t209 - t203;
t204 = t136 * t115;
t208 = t133 * t116;
t62 = t130 * t204 - t208;
t239 = -g(2) * t60 - g(3) * t62 + t115 * t232;
t238 = t54 ^ 2;
t125 = t130 ^ 2;
t97 = -qJD(4) + t105;
t5 = -t97 * pkin(4) - t242;
t236 = t5 + t242;
t233 = pkin(7) * t129;
t231 = g(2) * t133;
t230 = g(3) * t136;
t149 = qJD(1) * t79;
t51 = t129 * t149;
t68 = pkin(3) * t177 + qJ(2) * t196;
t36 = t51 * pkin(4) + qJD(5) + t68;
t229 = t36 * t54;
t228 = t54 * t51;
t227 = t134 * t38 - t31;
t76 = t135 * t85;
t42 = -t180 + t76 + (-pkin(3) - t219) * t130;
t213 = t129 * t132;
t220 = t132 * t85 + t104;
t47 = -pkin(7) * t213 + t220;
t226 = t131 * t42 + t134 * t47;
t210 = t131 * t132;
t78 = t205 - t210;
t225 = (t181 - t190) * t78;
t224 = t130 * t149 - t244;
t223 = t181 * t158;
t33 = t134 * t39;
t222 = t51 * qJ(5);
t193 = qJD(3) * t135;
t194 = qJD(2) * t130;
t221 = t135 * t194 + t85 * t193;
t217 = qJD(3) * t67;
t215 = (-qJ(5) - pkin(7) - pkin(6)) * t129;
t137 = qJD(1) ^ 2;
t214 = t124 * t137;
t211 = t130 * t133;
t207 = t133 * t132;
t206 = t133 * t135;
t202 = t136 * t132;
t201 = t136 * t135;
t84 = t135 * pkin(3) + pkin(4) * t116;
t100 = t129 * pkin(3) * t193;
t74 = t129 * qJD(2) + t100;
t106 = pkin(3) * t213;
t77 = t129 * qJ(2) + t106;
t200 = t136 * pkin(1) + t133 * qJ(2);
t198 = t124 + t125;
t127 = t135 ^ 2;
t197 = t132 ^ 2 - t127;
t192 = qJD(4) * t131;
t191 = qJD(4) * t134;
t189 = qJD(3) + t105;
t187 = qJD(1) * qJD(3);
t184 = qJDD(1) * t135;
t182 = t130 * qJDD(1);
t178 = qJ(2) * qJD(3) * t130;
t175 = qJ(2) * t182;
t174 = t135 * t188;
t173 = t135 * t187;
t171 = t131 * t184;
t103 = -qJDD(3) + t182;
t157 = qJD(1) * t178;
t66 = t85 * qJDD(1) + qJDD(2);
t58 = t135 * t66;
t12 = -t103 * pkin(3) + t58 + (-pkin(7) * t183 - t157) * t135 + (-t175 - t217 + (qJD(3) * t233 - t194) * qJD(1)) * t132;
t164 = t130 * t174 + t132 * t66 + t135 * t175 + t67 * t193;
t143 = -t132 * t157 + t164;
t18 = (-t173 - t185) * t233 + t143;
t170 = t134 * t12 - t131 * t18;
t34 = t148 * qJD(3) + t221;
t35 = -t132 * t194 + (-t104 + (-t85 + t233) * t132) * qJD(3);
t168 = -t131 * t34 + t134 * t35;
t167 = -t131 * t38 - t33;
t166 = -t131 * t47 + t134 * t42;
t165 = -t131 * t12 - t134 * t18 - t24 * t191 + t39 * t192;
t163 = t198 * t137;
t44 = qJ(2) * t183 + qJD(1) * t100 + qJDD(1) * t106 + t129 * t188;
t162 = qJD(1) * t189;
t161 = t181 * t135;
t160 = t103 + t182;
t159 = 0.2e1 * t198;
t155 = -t230 + t231;
t154 = -t131 * t24 - t33;
t152 = qJD(3) * (t105 + t190);
t150 = t131 * t35 + t134 * t34 + t42 * t191 - t47 * t192;
t20 = (t171 + (qJD(1) * t161 + t185) * t134) * t129 - t223;
t147 = t20 * pkin(4) + qJDD(5) + t44;
t146 = t216 + t243;
t145 = -t105 ^ 2 - t214;
t144 = t159 * t188 + t230;
t142 = t154 * qJD(4) + t170;
t61 = t130 * t208 - t204;
t63 = t130 * t203 + t209;
t141 = g(2) * t61 - g(3) * t63 + t116 * t232 + t68 * t51 + t165;
t139 = -t68 * t54 + t142 + t239;
t118 = t133 * pkin(1);
t111 = t134 * pkin(3) + pkin(4);
t91 = -qJDD(4) + t103;
t83 = t132 * pkin(3) + pkin(4) * t115;
t80 = pkin(2) + t84;
t72 = t130 * t201 + t207;
t71 = t130 * t202 - t206;
t70 = t130 * t206 - t202;
t69 = -t130 * t207 - t201;
t65 = t78 * t129;
t64 = t79 * t129;
t49 = t51 ^ 2;
t30 = (t134 * t161 - t181 * t210) * t129;
t29 = t244 * t129;
t21 = -t49 + t238;
t17 = -qJ(5) * t64 + t226;
t15 = -t130 * pkin(4) - t65 * qJ(5) + t166;
t14 = -t54 * t97 + (-t171 + (-t181 * t195 - t185) * t134) * t129 + t223;
t13 = -t51 * t97 - t19;
t9 = -t48 + t227;
t8 = t167 + t222;
t7 = -t154 - t222;
t4 = t29 * qJ(5) - qJD(4) * t226 - t65 * qJD(5) + t168;
t3 = -t30 * qJ(5) - t64 * qJD(5) + t150;
t2 = -qJ(5) * t20 - qJD(5) * t51 - t165;
t1 = -t91 * pkin(4) + t19 * qJ(5) - t54 * qJD(5) + t142;
t6 = [qJDD(1), t240, t155, t146 * t130, -t146 * t129, t159 * t186 + t144 - t231, -t112 * pkin(1) - g(2) * t200 - g(3) * t118 + (t198 * t186 + t144) * qJ(2), (qJDD(1) * t127 - 0.2e1 * t132 * t173) * t124, (-t132 * t184 + t187 * t197) * t237, (t132 * t152 - t135 * t160) * t129, (t132 * t160 + t135 * t152) * t129, t103 * t130, -g(2) * t72 - g(3) * t70 - t76 * t103 - t58 * t130 + (t105 * t130 + (t237 + t125) * qJD(1)) * qJ(2) * t193 + (qJD(3) * t85 * t105 + t151 * t237 + (qJ(2) * t103 + qJD(2) * t105 + t217 + t241) * t130) * t132, (-t132 * t178 + t221) * t105 + t220 * t103 + t143 * t130 + g(2) * t71 - g(3) * t69 + (t174 + (-t132 * t187 + t184) * qJ(2)) * t237, -t19 * t65 - t29 * t54, t19 * t64 - t20 * t65 + t29 * t51 - t30 * t54, t19 * t130 + t29 * t97 - t65 * t91, t20 * t130 + t30 * t97 + t64 * t91, t91 * t130, -t168 * t97 - t166 * t91 - t170 * t130 + t74 * t51 + t77 * t20 + t44 * t64 + t68 * t30 - g(2) * t63 - g(3) * t61 + (-t130 * t154 + t226 * t97) * qJD(4), g(2) * t62 - g(3) * t60 - t165 * t130 + t150 * t97 - t77 * t19 + t226 * t91 - t68 * t29 + t44 * t65 + t74 * t54, -t1 * t65 + t129 * t240 + t15 * t19 - t17 * t20 - t2 * t64 + t5 * t29 - t3 * t51 - t7 * t30 - t4 * t54, t2 * t17 + t7 * t3 + t1 * t15 + t5 * t4 + t147 * (t64 * pkin(4) + t77) + t36 * (t30 * pkin(4) + t74) - g(2) * (t133 * t83 + t200) - g(3) * (-t133 * t215 + t211 * t80 + t118) + (-g(2) * (t130 * t80 - t215) - g(3) * (-qJ(2) - t83)) * t136; 0, 0, 0, -t182, t183, -t163, -qJ(2) * t163 - t243, 0, 0, 0, 0, 0, -t135 * t103 + t132 * t145, t132 * t103 + t135 * t145, 0, 0, 0, 0, 0, -t51 * t196 - t224 * t97 - t78 * t91, -t54 * t196 + t225 * t97 + t79 * t91, t19 * t78 - t20 * t79 - t224 * t54 - t225 * t51, t1 * t78 - t36 * t196 + t2 * t79 + t224 * t5 + t225 * t7 - t240; 0, 0, 0, 0, 0, 0, 0, t135 * t132 * t214, -t197 * t214, (-t132 * t162 + t184) * t129, (-t135 * t162 - t185) * t129, -t103, -g(2) * t69 - g(3) * t71 + t58 + (-t130 * t162 - t214) * t218 + (-t189 * t67 + t232 - t241) * t132, g(1) * t212 + g(2) * t70 - g(3) * t72 - t59 * t105 + (t189 * t190 + t214) * t219 - t164, t228, t21, t13, t14, -t91, t167 * t97 + (-t134 * t91 - t176 * t51 + t192 * t97) * pkin(3) + t139, -t227 * t97 + (t131 * t91 - t176 * t54 + t191 * t97) * pkin(3) + t141, t111 * t19 + (t7 + t8) * t54 + (-t5 + t9) * t51 + (-t131 * t20 + (t131 * t54 - t134 * t51) * qJD(4)) * pkin(3), t1 * t111 - t7 * t9 - t5 * t8 - pkin(4) * t229 + t83 * t232 - g(2) * (-t136 * t84 - t211 * t83) - g(3) * (t136 * t130 * t83 - t133 * t84) + (-t36 * t176 + t2 * t131 + (-t131 * t5 + t134 * t7) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t21, t13, t14, -t91, t154 * t97 + t139, -t169 * t97 + t141, pkin(4) * t19 - t236 * t51, t236 * t7 + (t1 - t229 + t239) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49 - t238, g(1) * t130 - t129 * t155 + t5 * t54 + t7 * t51 + t147;];
tau_reg = t6;
