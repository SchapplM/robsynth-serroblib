% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR8
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:16
% EndTime: 2019-12-31 20:18:23
% DurationCPUTime: 2.80s
% Computational Cost: add. (3917->306), mult. (9486->412), div. (0->0), fcn. (7410->12), ass. (0->178)
t164 = cos(qJ(5));
t210 = qJD(5) * t164;
t157 = sin(pkin(9));
t158 = cos(pkin(9));
t162 = sin(qJ(2));
t166 = cos(qJ(2));
t121 = t157 * t166 + t158 * t162;
t112 = t121 * qJD(1);
t161 = sin(qJ(4));
t120 = -t157 * t162 + t158 * t166;
t110 = t120 * qJD(1);
t165 = cos(qJ(4));
t96 = t165 * t110;
t68 = -t112 * t161 + t96;
t262 = t164 * t68;
t270 = t210 - t262;
t152 = qJ(2) + pkin(9) + qJ(4);
t146 = sin(t152);
t163 = sin(qJ(1));
t167 = cos(qJ(1));
t188 = g(1) * t167 + g(2) * t163;
t269 = t188 * t146;
t154 = qJD(2) + qJD(4);
t225 = t154 * t68;
t212 = qJD(4) * t161;
t111 = t121 * qJD(2);
t77 = -qJD(1) * t111 + t120 * qJDD(1);
t209 = qJD(1) * qJD(2);
t201 = t166 * t209;
t202 = t162 * t209;
t78 = t121 * qJDD(1) - t157 * t202 + t158 * t201;
t25 = qJD(4) * t96 - t112 * t212 + t161 * t77 + t165 * t78;
t268 = t25 - t225;
t153 = qJDD(2) + qJDD(4);
t160 = sin(qJ(5));
t185 = t110 * t161 + t165 * t112;
t211 = qJD(5) * t160;
t13 = t160 * t153 + t154 * t210 + t164 * t25 - t185 * t211;
t54 = t154 * t160 + t164 * t185;
t14 = qJD(5) * t54 - t164 * t153 + t160 * t25;
t52 = -t164 * t154 + t160 * t185;
t267 = t13 * t164 - t160 * t14 - t270 * t52;
t11 = t13 * t160;
t266 = t270 * t54 + t11;
t26 = t185 * qJD(4) + t161 * t78 - t165 * t77;
t24 = qJDD(5) + t26;
t18 = t160 * t24;
t233 = t54 * t185;
t260 = qJD(5) - t68;
t57 = t260 * t210;
t265 = -t260 * t262 + t18 - t233 + t57;
t241 = pkin(7) * t112;
t232 = qJ(3) + pkin(6);
t136 = t232 * t166;
t128 = qJD(1) * t136;
t115 = t157 * t128;
t135 = t232 * t162;
t127 = qJD(1) * t135;
t227 = qJD(2) * pkin(2);
t119 = -t127 + t227;
t75 = t158 * t119 - t115;
t47 = qJD(2) * pkin(3) - t241 + t75;
t242 = pkin(7) * t110;
t220 = t158 * t128;
t76 = t157 * t119 + t220;
t51 = t76 + t242;
t28 = -t161 * t51 + t165 * t47;
t20 = -pkin(4) * t154 - t28;
t264 = t20 * t68;
t147 = cos(t152);
t238 = g(3) * t147;
t193 = qJD(2) * t232;
t108 = -qJD(3) * t162 - t166 * t193;
t74 = qJDD(2) * pkin(2) + t108 * qJD(1) - qJDD(1) * t135;
t107 = qJD(3) * t166 - t162 * t193;
t81 = t107 * qJD(1) + qJDD(1) * t136;
t40 = -t157 * t81 + t158 * t74;
t27 = qJDD(2) * pkin(3) - pkin(7) * t78 + t40;
t29 = t161 * t47 + t165 * t51;
t41 = t157 * t74 + t158 * t81;
t30 = pkin(7) * t77 + t41;
t246 = t29 * qJD(4) + t161 * t30 - t165 * t27;
t3 = -pkin(4) * t153 + t246;
t263 = t3 + t238;
t261 = t185 * t68;
t226 = t154 * t185;
t258 = -t26 + t226;
t256 = t185 ^ 2 - t68 ^ 2;
t140 = g(3) * t146;
t245 = (qJD(4) * t47 + t30) * t165 + t161 * t27 - t51 * t212;
t149 = pkin(2) * t166 + pkin(1);
t129 = -t149 * qJD(1) + qJD(3);
t84 = -pkin(3) * t110 + t129;
t255 = t188 * t147 - t84 * t68 + t140 - t245;
t253 = pkin(4) * t185;
t148 = pkin(2) * t158 + pkin(3);
t243 = pkin(2) * t157;
t214 = t161 * t148 + t165 * t243;
t106 = pkin(8) + t214;
t87 = t162 * qJD(1) * pkin(2) + pkin(3) * t112;
t252 = (-pkin(8) * t68 + qJD(5) * t106 + t253 + t87) * t260;
t251 = (t260 * pkin(8) + t253) * t260;
t234 = t52 * t185;
t250 = t260 * t185;
t21 = pkin(8) * t154 + t29;
t31 = -pkin(4) * t68 - pkin(8) * t185 + t84;
t6 = -t160 * t21 + t164 * t31;
t249 = t164 * t269 - t6 * t185 + t20 * t211;
t7 = t160 * t31 + t164 * t21;
t248 = t263 * t160 + t7 * t185 + t20 * t210;
t247 = -t185 * t84 - t238 - t246 + t269;
t184 = t165 * t120 - t121 * t161;
t199 = pkin(8) * t153 + qJD(5) * t31 + t245;
t80 = t120 * t161 + t121 * t165;
t90 = -pkin(3) * t120 - t149;
t35 = -pkin(4) * t184 - pkin(8) * t80 + t90;
t85 = -t158 * t135 - t136 * t157;
t60 = -pkin(7) * t121 + t85;
t86 = -t157 * t135 + t158 * t136;
t61 = pkin(7) * t120 + t86;
t37 = t161 * t60 + t165 * t61;
t114 = t120 * qJD(2);
t42 = t184 * qJD(4) - t111 * t161 + t114 * t165;
t36 = t161 * t61 - t165 * t60;
t58 = -t107 * t157 + t158 * t108;
t44 = -pkin(7) * t114 + t58;
t59 = t158 * t107 + t157 * t108;
t45 = -pkin(7) * t111 + t59;
t8 = -t36 * qJD(4) + t161 * t44 + t165 * t45;
t244 = -(qJD(5) * t35 + t8) * t260 + t199 * t184 + t20 * t42 - t37 * t24 + t3 * t80;
t237 = g(3) * t166;
t236 = t20 * t80;
t235 = t35 * t24;
t181 = t148 * t165 - t161 * t243;
t82 = t127 * t157 - t220;
t55 = t82 - t242;
t83 = -t158 * t127 - t115;
t56 = t83 - t241;
t229 = -t181 * qJD(4) + t161 * t55 + t165 * t56;
t228 = t214 * qJD(4) - t161 * t56 + t165 * t55;
t224 = t160 * t54;
t219 = t160 * t163;
t218 = t160 * t167;
t217 = t163 * t164;
t216 = t164 * t167;
t155 = t162 ^ 2;
t213 = -t166 ^ 2 + t155;
t208 = t166 * qJDD(1);
t151 = t162 * t227;
t204 = t80 * t211;
t88 = pkin(3) * t111 + t151;
t176 = pkin(2) * t202 - t149 * qJDD(1) + qJDD(3);
t50 = -pkin(3) * t77 + t176;
t5 = pkin(4) * t26 - pkin(8) * t25 + t50;
t198 = qJD(5) * t21 - t5;
t190 = t160 * t260;
t187 = g(1) * t163 - g(2) * t167;
t186 = t24 * t80 + t260 * t42;
t19 = t164 * t24;
t183 = t19 - (-t160 * t68 + t211) * t260;
t182 = -t199 + t140;
t179 = -0.2e1 * pkin(1) * t209 - pkin(6) * qJDD(2);
t178 = -pkin(8) * t24 + t260 * t28 - t264;
t175 = -t106 * t24 + t229 * t260 - t264;
t168 = qJD(2) ^ 2;
t174 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t168 + t187;
t169 = qJD(1) ^ 2;
t173 = pkin(1) * t169 - pkin(6) * qJDD(1) + t188;
t105 = -pkin(4) - t181;
t104 = t147 * t216 + t219;
t103 = -t147 * t218 + t217;
t102 = -t147 * t217 + t218;
t101 = t147 * t219 + t216;
t43 = t80 * qJD(4) + t165 * t111 + t114 * t161;
t15 = pkin(4) * t43 - pkin(8) * t42 + t88;
t9 = t37 * qJD(4) + t161 * t45 - t165 * t44;
t4 = t164 * t5;
t1 = [qJDD(1), t187, t188, qJDD(1) * t155 + 0.2e1 * t162 * t201, 0.2e1 * t162 * t208 - 0.2e1 * t213 * t209, qJDD(2) * t162 + t166 * t168, qJDD(2) * t166 - t162 * t168, 0, t179 * t162 + t174 * t166, -t174 * t162 + t179 * t166, t110 * t59 - t111 * t76 - t112 * t58 - t114 * t75 + t120 * t41 - t121 * t40 + t77 * t86 - t78 * t85 - t188, t41 * t86 + t76 * t59 + t40 * t85 + t75 * t58 - t176 * t149 + t129 * t151 - g(1) * (-t149 * t163 + t167 * t232) - g(2) * (t149 * t167 + t163 * t232), t185 * t42 + t25 * t80, t184 * t25 - t185 * t43 - t26 * t80 + t42 * t68, t153 * t80 + t154 * t42, t153 * t184 - t154 * t43, 0, t147 * t187 - t153 * t36 - t154 * t9 - t184 * t50 + t26 * t90 + t43 * t84 - t68 * t88, -t146 * t187 - t153 * t37 - t154 * t8 + t185 * t88 + t25 * t90 + t42 * t84 + t50 * t80, -t54 * t204 + (t13 * t80 + t42 * t54) * t164, (-t164 * t52 - t224) * t42 + (-t11 - t14 * t164 + (t160 * t52 - t164 * t54) * qJD(5)) * t80, -t13 * t184 + t164 * t186 - t204 * t260 + t43 * t54, t14 * t184 - t160 * t186 - t43 * t52 - t57 * t80, -t184 * t24 + t260 * t43, -g(1) * t102 - g(2) * t104 + t36 * t14 - t4 * t184 + t6 * t43 + t9 * t52 + (t15 * t260 + t235 + (t184 * t21 - t260 * t37 + t236) * qJD(5)) * t164 + t244 * t160, -g(1) * t101 - g(2) * t103 + t36 * t13 - t7 * t43 + t9 * t54 + (-(-qJD(5) * t37 + t15) * t260 - t235 - t198 * t184 - qJD(5) * t236) * t160 + t244 * t164; 0, 0, 0, -t162 * t169 * t166, t213 * t169, t162 * qJDD(1), t208, qJDD(2), t173 * t162 - t237, g(3) * t162 + t173 * t166, (t76 + t82) * t112 + (t75 - t83) * t110 + (t157 * t77 - t158 * t78) * pkin(2), -t75 * t82 - t76 * t83 + (-t237 + t157 * t41 + t158 * t40 + (-qJD(1) * t129 + t188) * t162) * pkin(2), -t261, t256, t268, t258, t153, t153 * t181 - t154 * t228 + t68 * t87 + t247, -t153 * t214 + t154 * t229 - t185 * t87 + t255, t266, -t224 * t260 + t267, t265, t183 + t234, -t250, t105 * t14 + t228 * t52 + (-t263 - t252) * t164 + t175 * t160 + t249, t105 * t13 + t228 * t54 + t175 * t164 + (-t269 + t252) * t160 + t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110 ^ 2 - t112 ^ 2, -t110 * t76 + t112 * t75 + t176 - t187, 0, 0, 0, 0, 0, t26 + t226, t25 + t225, 0, 0, 0, 0, 0, t183 - t234, -t164 * t260 ^ 2 - t18 - t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261, t256, t268, t258, t153, t154 * t29 + t247, t154 * t28 + t255, t266, -t190 * t54 + t267, t265, -t190 * t260 + t19 + t234, -t250, -pkin(4) * t14 - t29 * t52 + t178 * t160 + (-t263 - t251) * t164 + t249, -pkin(4) * t13 - t29 * t54 + t178 * t164 + (-t269 + t251) * t160 + t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t52, -t52 ^ 2 + t54 ^ 2, t260 * t52 + t13, t260 * t54 - t14, t24, -g(1) * t103 + g(2) * t101 + t160 * t182 - t20 * t54 - t21 * t210 + t260 * t7 + t4, g(1) * t104 - g(2) * t102 + t160 * t198 + t164 * t182 + t20 * t52 + t260 * t6;];
tau_reg = t1;
