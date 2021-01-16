% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:18
% EndTime: 2021-01-15 20:18:28
% DurationCPUTime: 2.01s
% Computational Cost: add. (3362->312), mult. (7866->387), div. (0->0), fcn. (5721->12), ass. (0->170)
t155 = sin(qJ(4));
t153 = cos(pkin(8));
t152 = sin(pkin(8));
t158 = cos(qJ(2));
t204 = qJD(1) * qJD(2);
t196 = t158 * t204;
t156 = sin(qJ(2));
t197 = t156 * t204;
t203 = t156 * qJDD(1);
t199 = -t153 * t197 + (-t196 - t203) * t152;
t202 = t158 * qJDD(1);
t175 = t153 * t202 + t199;
t231 = cos(qJ(4));
t211 = t153 * t158;
t185 = t152 * t156 - t211;
t94 = t185 * qJD(1);
t103 = t152 * t158 + t153 * t156;
t96 = t103 * qJD(1);
t181 = t155 * t94 - t231 * t96;
t171 = qJDD(1) * t103 - t152 * t197;
t62 = t153 * t196 + t171;
t165 = qJD(4) * t181 - t155 * t62 + t231 * t175;
t148 = qJD(2) + qJD(4);
t217 = t181 * t148;
t247 = t165 - t217;
t198 = qJD(4) * t231;
t205 = qJD(4) * t155;
t179 = t155 * t175 - t94 * t198 - t96 * t205 + t231 * t62;
t52 = -t155 * t96 - t231 * t94;
t218 = t148 * t52;
t7 = t179 - t218;
t224 = t52 ^ 2;
t245 = t181 ^ 2;
t246 = -t224 + t245;
t226 = t158 * pkin(2);
t136 = pkin(1) + t226;
t112 = -t136 * qJD(1) + qJD(3);
t68 = pkin(3) * t94 + t112;
t20 = -pkin(4) * t52 + qJ(5) * t181 + t68;
t244 = t20 * t52;
t243 = t68 * t52;
t223 = t181 * t52;
t242 = t181 * t68;
t149 = qJ(2) + pkin(8);
t143 = cos(t149);
t144 = qJ(4) + t149;
t133 = sin(t144);
t134 = cos(t144);
t209 = t134 * pkin(4) + t133 * qJ(5);
t241 = pkin(3) * t143 + t209 + t226;
t146 = qJDD(2) + qJDD(4);
t141 = t146 * pkin(4);
t236 = qJDD(5) - t141;
t240 = -t181 * t20 + t236;
t27 = -pkin(4) * t181 - qJ(5) * t52;
t135 = pkin(2) * t153 + pkin(3);
t227 = t152 * pkin(2);
t200 = t155 * t227;
t233 = pkin(7) * t94;
t154 = -qJ(3) - pkin(6);
t120 = t154 * t156;
t108 = qJD(1) * t120;
t121 = t154 * t158;
t109 = qJD(1) * t121;
t212 = t153 * t109;
t66 = -t108 * t152 + t212;
t41 = t66 + t233;
t232 = pkin(7) * t96;
t98 = t152 * t109;
t67 = t153 * t108 + t98;
t42 = t67 - t232;
t239 = qJD(4) * t200 - t135 * t198 + t155 * t41 + t231 * t42;
t208 = t155 * t135 + t231 * t227;
t137 = t146 * qJ(5);
t140 = t148 * qJD(5);
t238 = t137 + t140;
t157 = sin(qJ(1));
t159 = cos(qJ(1));
t237 = g(1) * t157 - g(2) * t159;
t69 = t153 * t120 + t121 * t152;
t45 = -pkin(7) * t103 + t69;
t70 = t152 * t120 - t153 * t121;
t46 = -pkin(7) * t185 + t70;
t26 = t155 * t45 + t231 * t46;
t182 = -t155 * t46 + t231 * t45;
t176 = t185 * qJD(2);
t191 = qJD(2) * t154;
t91 = qJD(3) * t158 + t156 * t191;
t92 = -qJD(3) * t156 + t158 * t191;
t43 = -t152 * t91 + t153 * t92;
t34 = pkin(7) * t176 + t43;
t177 = t103 * qJD(2);
t44 = t152 * t92 + t153 * t91;
t35 = -pkin(7) * t177 + t44;
t4 = t182 * qJD(4) + t155 * t34 + t231 * t35;
t235 = t133 * t237 + t146 * t26 + t148 * t4;
t230 = pkin(2) * t156;
t228 = g(3) * t158;
t59 = qJDD(2) * pkin(2) + t92 * qJD(1) + qJDD(1) * t120;
t65 = t91 * qJD(1) - qJDD(1) * t121;
t29 = t152 * t59 + t153 * t65;
t222 = qJD(5) - t239;
t221 = t208 * qJD(4) - t155 * t42 + t231 * t41;
t220 = qJD(2) * pkin(2);
t102 = t108 + t220;
t60 = t153 * t102 + t98;
t38 = qJD(2) * pkin(3) - t232 + t60;
t61 = t152 * t102 - t212;
t40 = t61 - t233;
t18 = t155 * t38 + t231 * t40;
t219 = t148 * t18;
t216 = t133 * t157;
t215 = t133 * t159;
t214 = t134 * t157;
t213 = t134 * t159;
t17 = -t155 * t40 + t231 * t38;
t210 = qJD(5) - t17;
t150 = t156 ^ 2;
t207 = -t158 ^ 2 + t150;
t206 = qJD(1) * t156;
t201 = pkin(2) * t197 + qJDD(3);
t139 = t156 * t220;
t71 = pkin(2) * t206 + pkin(3) * t96;
t142 = sin(t149);
t194 = -pkin(3) * t142 - pkin(4) * t133 - t230;
t28 = -t152 * t65 + t153 * t59;
t16 = qJDD(2) * pkin(3) - pkin(7) * t62 + t28;
t19 = pkin(7) * t175 + t29;
t190 = t155 * t16 + t231 * t19 + t38 * t198 - t40 * t205;
t189 = t155 * t19 - t231 * t16 + t40 * t198 + t38 * t205;
t188 = g(1) * t159 + g(2) * t157;
t183 = pkin(1) + t241;
t180 = -0.2e1 * pkin(1) * t204 - pkin(6) * qJDD(2);
t178 = g(1) * t213 + g(2) * t214 + g(3) * t133 - t190;
t174 = t231 * t135 - t200;
t88 = -t136 * qJDD(1) + t201;
t173 = g(1) * t215 + g(2) * t216 - g(3) * t134 - t189;
t172 = t231 * t185;
t170 = t179 + t218;
t169 = t148 * t17 + t178;
t5 = t26 * qJD(4) + t155 * t35 - t231 * t34;
t168 = g(1) * t214 - g(2) * t213 + t146 * t182 - t148 * t5;
t72 = pkin(3) * t177 + t139;
t160 = qJD(2) ^ 2;
t167 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t160 + t237;
t161 = qJD(1) ^ 2;
t166 = pkin(1) * t161 - pkin(6) * qJDD(1) + t188;
t75 = pkin(3) * t185 - t136;
t64 = t231 * t103 - t155 * t185;
t164 = -t221 * t148 + t173;
t163 = -t173 + t240;
t162 = -t165 - t217;
t39 = -t199 * pkin(3) + (-pkin(1) + (-pkin(3) * t153 - pkin(2)) * t158) * qJDD(1) + t201;
t3 = -pkin(4) * t165 - qJ(5) * t179 + qJD(5) * t181 + t39;
t147 = -pkin(7) + t154;
t114 = qJ(5) * t213;
t113 = qJ(5) * t214;
t89 = -pkin(4) - t174;
t87 = qJ(5) + t208;
t63 = t103 * t155 + t172;
t31 = t64 * qJD(4) - t155 * t176 + t231 * t177;
t30 = t103 * t205 + t148 * t172 + t155 * t177;
t24 = t63 * pkin(4) - t64 * qJ(5) + t75;
t23 = t27 + t71;
t9 = t148 * qJ(5) + t18;
t8 = -t148 * pkin(4) + t210;
t6 = t31 * pkin(4) + t30 * qJ(5) - t64 * qJD(5) + t72;
t2 = t189 + t236;
t1 = t190 + t238;
t10 = [qJDD(1), t237, t188, qJDD(1) * t150 + 0.2e1 * t156 * t196, 0.2e1 * t156 * t202 - 0.2e1 * t207 * t204, qJDD(2) * t156 + t158 * t160, qJDD(2) * t158 - t156 * t160, 0, t156 * t180 + t158 * t167, -t156 * t167 + t158 * t180, t136 * t175 + t88 * t185 + t69 * qJDD(2) + t237 * t143 + (t103 * t112 + t94 * t230 + t43) * qJD(2), -t70 * qJDD(2) + t88 * t103 - t136 * t62 - t237 * t142 + (-t112 * t185 + t96 * t230 - t44) * qJD(2), -t44 * t94 + t70 * t175 - t29 * t185 - t43 * t96 - t69 * t62 - t28 * t103 + (-t103 * t61 + t185 * t60) * qJD(2) - t188, t29 * t70 + t61 * t44 + t28 * t69 + t60 * t43 - t88 * t136 + t112 * t139 - g(1) * (-t136 * t157 - t154 * t159) - g(2) * (t136 * t159 - t154 * t157), t179 * t64 + t181 * t30, t165 * t64 - t179 * t63 + t181 * t31 - t30 * t52, t146 * t64 - t148 * t30, -t146 * t63 - t148 * t31, 0, -t165 * t75 + t31 * t68 + t39 * t63 - t52 * t72 + t168, t179 * t75 - t181 * t72 - t30 * t68 + t39 * t64 - t235, -t165 * t24 + t20 * t31 + t3 * t63 - t52 * t6 + t168, -t1 * t63 + t165 * t26 - t179 * t182 - t181 * t5 + t2 * t64 - t30 * t8 - t31 * t9 + t4 * t52 - t188, -t179 * t24 + t181 * t6 + t20 * t30 - t3 * t64 + t235, t1 * t26 - t2 * t182 + t20 * t6 + t3 * t24 + t9 * t4 + t8 * t5 + (g(1) * t147 - g(2) * t183) * t159 + (g(1) * t183 + g(2) * t147) * t157; 0, 0, 0, -t156 * t161 * t158, t207 * t161, t203, t202, qJDD(2), t156 * t166 - t228, g(3) * t156 + t158 * t166, -g(3) * t143 - qJD(2) * t66 - t112 * t96 + t188 * t142 + (qJDD(2) * t153 - t94 * t206) * pkin(2) + t28, g(3) * t142 + qJD(2) * t67 + t112 * t94 + t188 * t143 + (-qJDD(2) * t152 - t96 * t206) * pkin(2) - t29, (t61 + t66) * t96 + (t67 - t60) * t94 + (t152 * t175 - t153 * t62) * pkin(2), -t60 * t66 - t61 * t67 + (-t228 + t29 * t152 + t153 * t28 + (-qJD(1) * t112 + t188) * t156) * pkin(2), t223, t246, t7, t247, t146, t146 * t174 + t52 * t71 + t164 + t242, -t208 * t146 + t239 * t148 + t181 * t71 + t178 - t243, -t146 * t89 + t23 * t52 + t164 - t240, t165 * t87 + t179 * t89 + (t222 - t8) * t52 + (-t221 - t9) * t181, t146 * t87 + t222 * t148 - t181 * t23 - t178 + t238 + t244, t1 * t87 + t2 * t89 - t20 * t23 - g(1) * (t159 * t194 + t114) - g(2) * (t157 * t194 + t113) - g(3) * t241 + t222 * t9 + t221 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t96 - t175, (qJD(1) * t211 - t94) * qJD(2) + t171, -t94 ^ 2 - t96 ^ 2, t60 * t96 + t61 * t94 - t237 + t88, 0, 0, 0, 0, 0, t162, t170, t162, -t224 - t245, -t170, t181 * t8 - t9 * t52 - t237 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, t246, t7, t247, t146, t173 + t219 + t242, t169 - t243, t27 * t52 + t141 - t163 + t219, -pkin(4) * t179 + qJ(5) * t165 - (-t18 + t9) * t181 - (t8 - t210) * t52, -t181 * t27 + 0.2e1 * t137 + 0.2e1 * t140 - t169 + t244, t1 * qJ(5) - t2 * pkin(4) - t20 * t27 - t8 * t18 - g(1) * (-pkin(4) * t215 + t114) - g(2) * (-pkin(4) * t216 + t113) - g(3) * t209 + t210 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146 + t223, t7, -t148 ^ 2 - t245, -t148 * t9 + t163;];
tau_reg = t10;
