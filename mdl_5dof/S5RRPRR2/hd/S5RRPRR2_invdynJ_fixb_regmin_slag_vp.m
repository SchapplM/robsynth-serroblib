% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:30
% EndTime: 2019-12-05 18:28:36
% DurationCPUTime: 2.37s
% Computational Cost: add. (3715->277), mult. (9083->374), div. (0->0), fcn. (7183->14), ass. (0->164)
t168 = sin(qJ(5));
t172 = cos(qJ(5));
t165 = sin(pkin(9));
t166 = cos(pkin(9));
t170 = sin(qJ(2));
t174 = cos(qJ(2));
t126 = -t165 * t170 + t166 * t174;
t116 = t126 * qJD(1);
t127 = t165 * t174 + t166 * t170;
t118 = t127 * qJD(1);
t169 = sin(qJ(4));
t173 = cos(qJ(4));
t219 = qJD(4) * t173;
t220 = qJD(4) * t169;
t117 = t127 * qJD(2);
t87 = -qJD(1) * t117 + t126 * qJDD(1);
t216 = qJD(1) * qJD(2);
t212 = t174 * t216;
t213 = t170 * t216;
t88 = t127 * qJDD(1) - t165 * t213 + t166 * t212;
t23 = t116 * t219 - t118 * t220 + t169 * t87 + t173 * t88;
t193 = t116 * t169 + t173 * t118;
t24 = t193 * qJD(4) + t169 * t88 - t173 * t87;
t77 = t173 * t116 - t118 * t169;
t29 = t168 * t77 + t172 * t193;
t182 = -qJD(5) * t29 - t168 * t23 - t172 * t24;
t162 = qJD(2) + qJD(4);
t159 = qJD(5) + t162;
t229 = t159 * t29;
t248 = t182 + t229;
t33 = -t168 * t193 + t172 * t77;
t256 = t29 * t33;
t228 = t159 * t33;
t217 = qJD(5) * t172;
t218 = qJD(5) * t168;
t4 = -t168 * t24 + t172 * t23 - t193 * t218 + t217 * t77;
t252 = t4 - t228;
t253 = t29 ^ 2 - t33 ^ 2;
t160 = qJ(2) + pkin(9) + qJ(4);
t154 = qJ(5) + t160;
t148 = sin(t154);
t149 = cos(t154);
t161 = qJDD(2) + qJDD(4);
t237 = pkin(7) * t118;
t233 = qJ(3) + pkin(6);
t143 = t233 * t174;
t132 = qJD(1) * t143;
t121 = t165 * t132;
t142 = t233 * t170;
t131 = qJD(1) * t142;
t230 = qJD(2) * pkin(2);
t125 = -t131 + t230;
t85 = t166 * t125 - t121;
t51 = qJD(2) * pkin(3) - t237 + t85;
t238 = pkin(7) * t116;
t222 = t166 * t132;
t86 = t165 * t125 + t222;
t57 = t86 + t238;
t194 = -t169 * t51 - t173 * t57;
t204 = qJD(2) * t233;
t114 = -t170 * qJD(3) - t174 * t204;
t84 = qJDD(2) * pkin(2) + t114 * qJD(1) - qJDD(1) * t142;
t113 = t174 * qJD(3) - t170 * t204;
t91 = t113 * qJD(1) + qJDD(1) * t143;
t36 = -t165 * t91 + t166 * t84;
t25 = qJDD(2) * pkin(3) - pkin(7) * t88 + t36;
t37 = t165 * t84 + t166 * t91;
t26 = pkin(7) * t87 + t37;
t184 = t194 * qJD(4) - t169 * t26 + t173 * t25;
t2 = pkin(4) * t161 - pkin(8) * t23 + t184;
t171 = sin(qJ(1));
t175 = cos(qJ(1));
t200 = g(1) * t175 + g(2) * t171;
t240 = (qJD(4) * t51 + t26) * t173 + t169 * t25 - t57 * t220;
t3 = -pkin(8) * t24 + t240;
t155 = pkin(2) * t174 + pkin(1);
t137 = -t155 * qJD(1) + qJD(3);
t94 = -pkin(3) * t116 + t137;
t42 = -pkin(4) * t77 + t94;
t261 = -g(3) * t149 + t200 * t148 - t168 * t3 + t172 * t2 - t42 * t29;
t258 = pkin(8) * t77;
t14 = -t194 + t258;
t260 = g(3) * t148 + t14 * t218 + t200 * t149 - t42 * t33;
t226 = t162 * t77;
t259 = t23 - t226;
t257 = t193 * t77;
t227 = t162 * t193;
t255 = -t24 + t227;
t254 = t193 ^ 2 - t77 ^ 2;
t151 = sin(t160);
t152 = cos(t160);
t251 = g(3) * t151 + t200 * t152 - t94 * t77 - t240;
t250 = (-t14 * t159 - t2) * t168 + t260;
t208 = -t169 * t57 + t173 * t51;
t246 = pkin(8) * t193;
t13 = t208 - t246;
t11 = pkin(4) * t162 + t13;
t225 = t172 * t14;
t198 = -t168 * t11 - t225;
t249 = t198 * qJD(5) + t261;
t153 = pkin(2) * t166 + pkin(3);
t239 = pkin(2) * t165;
t201 = t173 * t153 - t169 * t239;
t92 = t131 * t165 - t222;
t58 = t92 - t238;
t93 = -t166 * t131 - t121;
t59 = t93 - t237;
t243 = -t201 * qJD(4) + t169 * t58 + t173 * t59;
t112 = t153 * t169 + t173 * t239;
t242 = -t112 * qJD(4) + t169 * t59 - t173 * t58;
t241 = -g(3) * t152 + t200 * t151 - t94 * t193 + t184;
t234 = g(3) * t174;
t95 = -t166 * t142 - t143 * t165;
t70 = -pkin(7) * t127 + t95;
t96 = -t165 * t142 + t166 * t143;
t71 = pkin(7) * t126 + t96;
t231 = t169 * t70 + t173 * t71;
t224 = t243 - t246;
t223 = t242 + t258;
t69 = t166 * t113 + t165 * t114;
t163 = t170 ^ 2;
t221 = -t174 ^ 2 + t163;
t215 = t174 * qJDD(1);
t157 = t170 * t230;
t98 = pkin(3) * t117 + t157;
t97 = t170 * qJD(1) * pkin(2) + pkin(3) * t118;
t211 = -qJD(5) * t11 - t3;
t206 = -t169 * t71 + t173 * t70;
t68 = -t113 * t165 + t166 * t114;
t199 = g(1) * t171 - g(2) * t175;
t90 = t126 * t169 + t127 * t173;
t17 = -pkin(8) * t90 + t206;
t192 = t173 * t126 - t127 * t169;
t18 = pkin(8) * t192 + t231;
t197 = -t168 * t18 + t17 * t172;
t196 = t168 * t17 + t172 * t18;
t40 = t168 * t90 - t172 * t192;
t41 = t168 * t192 + t172 * t90;
t100 = -pkin(3) * t126 - t155;
t191 = -0.2e1 * pkin(1) * t216 - pkin(6) * qJDD(2);
t120 = t126 * qJD(2);
t47 = -pkin(7) * t120 + t68;
t48 = -pkin(7) * t117 + t69;
t190 = t169 * t47 + t173 * t48 + t70 * t219 - t71 * t220;
t188 = pkin(2) * t213 - t155 * qJDD(1) + qJDD(3);
t176 = qJD(2) ^ 2;
t186 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t176 + t199;
t177 = qJD(1) ^ 2;
t185 = pkin(1) * t177 - pkin(6) * qJDD(1) + t200;
t53 = -pkin(3) * t87 + t188;
t183 = -t231 * qJD(4) - t169 * t48 + t173 * t47;
t158 = qJDD(5) + t161;
t111 = pkin(4) + t201;
t52 = -pkin(4) * t192 + t100;
t43 = pkin(4) * t193 + t97;
t39 = t90 * qJD(4) + t173 * t117 + t120 * t169;
t38 = t192 * qJD(4) - t117 * t169 + t120 * t173;
t27 = pkin(4) * t39 + t98;
t10 = pkin(4) * t24 + t53;
t9 = t41 * qJD(5) + t168 * t38 + t172 * t39;
t8 = -t40 * qJD(5) - t168 * t39 + t172 * t38;
t7 = -pkin(8) * t38 + t183;
t6 = -pkin(8) * t39 + t190;
t1 = [qJDD(1), t199, t200, qJDD(1) * t163 + 0.2e1 * t170 * t212, 0.2e1 * t170 * t215 - 0.2e1 * t221 * t216, qJDD(2) * t170 + t174 * t176, qJDD(2) * t174 - t170 * t176, 0, t191 * t170 + t186 * t174, -t186 * t170 + t191 * t174, t116 * t69 - t117 * t86 - t118 * t68 - t120 * t85 + t126 * t37 - t127 * t36 + t87 * t96 - t88 * t95 - t200, t37 * t96 + t86 * t69 + t36 * t95 + t85 * t68 - t188 * t155 + t137 * t157 - g(1) * (-t155 * t171 + t175 * t233) - g(2) * (t155 * t175 + t171 * t233), t193 * t38 + t23 * t90, t192 * t23 - t193 * t39 - t24 * t90 + t38 * t77, t161 * t90 + t162 * t38, t161 * t192 - t162 * t39, 0, t100 * t24 + t199 * t152 + t206 * t161 + t183 * t162 - t192 * t53 + t94 * t39 - t77 * t98, t100 * t23 - t199 * t151 - t231 * t161 - t190 * t162 + t193 * t98 + t94 * t38 + t53 * t90, t29 * t8 + t4 * t41, t182 * t41 - t29 * t9 + t33 * t8 - t4 * t40, t158 * t41 + t159 * t8, -t158 * t40 - t159 * t9, 0, -t27 * t33 - t52 * t182 + t10 * t40 + t42 * t9 + (-qJD(5) * t196 - t168 * t6 + t172 * t7) * t159 + t197 * t158 + t199 * t149, t27 * t29 + t52 * t4 + t10 * t41 + t42 * t8 - (qJD(5) * t197 + t168 * t7 + t172 * t6) * t159 - t196 * t158 - t199 * t148; 0, 0, 0, -t170 * t177 * t174, t221 * t177, t170 * qJDD(1), t215, qJDD(2), t185 * t170 - t234, g(3) * t170 + t185 * t174, (t86 + t92) * t118 + (t85 - t93) * t116 + (t165 * t87 - t166 * t88) * pkin(2), -t85 * t92 - t86 * t93 + (-t234 + t165 * t37 + t166 * t36 + (-qJD(1) * t137 + t200) * t170) * pkin(2), -t257, t254, t259, t255, t161, t201 * t161 + t242 * t162 + t77 * t97 + t241, -t112 * t161 + t243 * t162 - t97 * t193 + t251, -t256, t253, t252, t248, t158, (t111 * t172 - t112 * t168) * t158 + t43 * t33 + (t224 * t168 + t223 * t172) * t159 + ((-t111 * t168 - t112 * t172) * t159 + t198) * qJD(5) + t261, -t43 * t29 + (-t111 * t158 - t2 + (qJD(5) * t112 - t223) * t159) * t168 + (-t112 * t158 + (-qJD(5) * t111 + t224) * t159 + t211) * t172 + t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 ^ 2 - t118 ^ 2, -t116 * t86 + t118 * t85 + t188 - t199, 0, 0, 0, 0, 0, t24 + t227, t23 + t226, 0, 0, 0, 0, 0, -t182 + t229, t4 + t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t257, t254, t259, t255, t161, -t194 * t162 + t241, t162 * t208 + t251, -t256, t253, t252, t248, t158, -(-t13 * t168 - t225) * t159 + (t158 * t172 - t159 * t218 + t193 * t33) * pkin(4) + t249, (t13 * t159 + t211) * t172 + (-t158 * t168 - t159 * t217 - t193 * t29) * pkin(4) + t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256, t253, t252, t248, t158, -t159 * t198 + t249, (-t3 + (-qJD(5) + t159) * t11) * t172 + t250;];
tau_reg = t1;
