% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:42
% EndTime: 2019-03-09 02:18:50
% DurationCPUTime: 3.00s
% Computational Cost: add. (4759->324), mult. (10679->418), div. (0->0), fcn. (8582->18), ass. (0->193)
t163 = cos(qJ(6));
t220 = qJD(6) * t163;
t157 = cos(pkin(11));
t165 = cos(qJ(4));
t229 = t165 * t157;
t155 = sin(pkin(11));
t161 = sin(qJ(4));
t230 = t155 * t161;
t113 = -t229 + t230;
t103 = t113 * qJD(1);
t114 = t155 * t165 + t157 * t161;
t104 = t114 * qJD(1);
t160 = sin(qJ(5));
t164 = cos(qJ(5));
t62 = t164 * t103 + t104 * t160;
t278 = t163 * t62;
t286 = t220 + t278;
t152 = pkin(11) + qJ(4);
t148 = qJ(5) + t152;
t135 = sin(t148);
t154 = qJ(1) + pkin(10);
t145 = sin(t154);
t147 = cos(t154);
t193 = g(1) * t147 + g(2) * t145;
t285 = t193 * t135;
t153 = qJD(4) + qJD(5);
t243 = t153 * t62;
t222 = qJD(5) * t164;
t223 = qJD(5) * t160;
t209 = qJD(1) * t230;
t216 = t157 * qJDD(1);
t217 = t155 * qJDD(1);
t224 = qJD(4) * t165;
t225 = qJD(1) * t157;
t211 = t161 * t216 + t165 * t217 + t224 * t225;
t70 = -qJD(4) * t209 + t211;
t106 = t114 * qJD(4);
t188 = t161 * t217 - t165 * t216;
t71 = qJD(1) * t106 + t188;
t34 = -t103 * t222 - t104 * t223 - t160 * t71 + t164 * t70;
t284 = t34 + t243;
t149 = qJDD(4) + qJDD(5);
t159 = sin(qJ(6));
t182 = -t103 * t160 + t164 * t104;
t221 = qJD(6) * t159;
t14 = t159 * t149 + t153 * t220 + t163 * t34 - t182 * t221;
t58 = t153 * t159 + t163 * t182;
t15 = qJD(6) * t58 - t163 * t149 + t159 * t34;
t56 = -t163 * t153 + t159 * t182;
t283 = t14 * t163 - t159 * t15 - t286 * t56;
t12 = t14 * t159;
t282 = t286 * t58 + t12;
t249 = t58 * t182;
t276 = qJD(6) + t62;
t35 = t182 * qJD(5) + t160 * t70 + t164 * t71;
t32 = qJDD(6) + t35;
t29 = t159 * t32;
t59 = t276 * t220;
t281 = t276 * t278 - t249 + t29 + t59;
t156 = sin(pkin(10));
t133 = pkin(1) * t156 + qJ(3);
t122 = t133 * qJD(1);
t143 = t157 * qJD(2);
t81 = t143 + (-pkin(7) * qJD(1) - t122) * t155;
t92 = t155 * qJD(2) + t157 * t122;
t82 = pkin(7) * t225 + t92;
t183 = -t161 * t81 - t165 * t82;
t46 = -pkin(8) * t103 - t183;
t241 = t160 * t46;
t265 = -t161 * t82 + t165 * t81;
t45 = -pkin(8) * t104 + t265;
t44 = qJD(4) * pkin(4) + t45;
t24 = t164 * t44 - t241;
t20 = -pkin(5) * t153 - t24;
t280 = t20 * t62;
t136 = cos(t148);
t253 = g(3) * t136;
t111 = qJD(1) * qJD(3) + t133 * qJDD(1);
t141 = t157 * qJDD(2);
t77 = t141 + (-pkin(7) * qJDD(1) - t111) * t155;
t84 = t155 * qJDD(2) + t157 * t111;
t78 = pkin(7) * t216 + t84;
t201 = -t161 * t78 + t165 * t77;
t22 = qJDD(4) * pkin(4) - t70 * pkin(8) + t183 * qJD(4) + t201;
t184 = t161 * t77 + t165 * t78;
t23 = -t71 * pkin(8) + t265 * qJD(4) + t184;
t237 = t164 * t46;
t25 = t160 * t44 + t237;
t260 = t25 * qJD(5) + t160 * t23 - t164 * t22;
t3 = -t149 * pkin(5) + t260;
t279 = t3 + t253;
t277 = t182 * t62;
t236 = t182 * t153;
t274 = -t35 + t236;
t272 = t182 ^ 2 - t62 ^ 2;
t129 = g(3) * t135;
t258 = (qJD(5) * t44 + t23) * t164 + t160 * t22 - t46 * t223;
t158 = cos(pkin(10));
t137 = -pkin(1) * t158 - pkin(2);
t120 = -pkin(3) * t157 + t137;
t101 = qJD(1) * t120 + qJD(3);
t69 = t103 * pkin(4) + t101;
t271 = t193 * t136 + t69 * t62 + t129 - t258;
t269 = pkin(5) * t182;
t138 = pkin(4) * t160 + pkin(9);
t268 = (pkin(4) * t104 + pkin(9) * t62 + qJD(6) * t138 + t269) * t276;
t267 = (t276 * pkin(9) + t269) * t276;
t250 = t56 * t182;
t266 = t276 * t182;
t218 = qJDD(1) * t137;
t119 = qJDD(3) + t218;
t192 = g(1) * t145 - g(2) * t147;
t264 = t192 - t119;
t21 = pkin(9) * t153 + t25;
t36 = t62 * pkin(5) - pkin(9) * t182 + t69;
t6 = -t159 * t21 + t163 * t36;
t263 = t163 * t285 - t6 * t182 + t20 * t221;
t7 = t159 * t36 + t163 * t21;
t262 = t279 * t159 + t7 * t182 + t20 * t220;
t261 = -t69 * t182 - t253 - t260 + t285;
t105 = t113 * qJD(4);
t74 = t164 * t113 + t114 * t160;
t47 = -t74 * qJD(5) - t164 * t105 - t160 * t106;
t75 = -t113 * t160 + t114 * t164;
t189 = t276 * t47 + t32 * t75;
t212 = t75 * t221;
t259 = -t163 * t189 + t212 * t276;
t207 = t149 * pkin(9) + qJD(6) * t36 + t258;
t248 = pkin(7) + t133;
t107 = t248 * t155;
t108 = t248 * t157;
t199 = -t165 * t107 - t108 * t161;
t54 = -pkin(8) * t114 + t199;
t244 = -t161 * t107 + t165 * t108;
t55 = -pkin(8) * t113 + t244;
t38 = t160 * t54 + t164 * t55;
t80 = pkin(4) * t113 + t120;
t40 = pkin(5) * t74 - pkin(9) * t75 + t80;
t37 = t160 * t55 - t164 * t54;
t175 = qJD(3) * t229 - t107 * t224 + (-qJD(3) * t155 - qJD(4) * t108) * t161;
t49 = -t106 * pkin(8) + t175;
t169 = -t114 * qJD(3) - t244 * qJD(4);
t50 = t105 * pkin(8) + t169;
t8 = -t37 * qJD(5) + t160 * t50 + t164 * t49;
t257 = t20 * t47 - (qJD(6) * t40 + t8) * t276 - t207 * t74 + t3 * t75 - t38 * t32;
t256 = pkin(4) * t106;
t252 = t20 * t75;
t251 = t40 * t32;
t48 = t75 * qJD(5) - t160 * t105 + t164 * t106;
t247 = t14 * t74 + t58 * t48;
t242 = t159 * t58;
t234 = t145 * t159;
t233 = t145 * t163;
t232 = t147 * t159;
t231 = t147 * t163;
t226 = t155 ^ 2 + t157 ^ 2;
t99 = qJDD(1) * t120 + qJDD(3);
t53 = t71 * pkin(4) + t99;
t5 = t35 * pkin(5) - t34 * pkin(9) + t53;
t206 = qJD(6) * t21 - t5;
t197 = t159 * t276;
t26 = t160 * t45 + t237;
t194 = pkin(4) * t223 - t26;
t162 = sin(qJ(1));
t166 = cos(qJ(1));
t191 = g(1) * t162 - g(2) * t166;
t190 = -t74 * t15 - t48 * t56;
t187 = t149 * t75 + t153 * t47;
t83 = -t111 * t155 + t141;
t186 = -t83 * t155 + t84 * t157;
t185 = t155 * (-t122 * t155 + t143) - t157 * t92;
t30 = t163 * t32;
t181 = t30 - (t159 * t62 + t221) * t276;
t180 = -t207 + t129;
t178 = -pkin(9) * t32 + t24 * t276 + t280;
t176 = -t218 + t264;
t27 = t164 * t45 - t241;
t173 = -t138 * t32 + t280 + (-pkin(4) * t222 + t27) * t276;
t172 = -t189 * t159 - t75 * t59;
t146 = cos(t152);
t144 = sin(t152);
t139 = -pkin(4) * t164 - pkin(5);
t88 = t136 * t231 + t234;
t87 = -t136 * t232 + t233;
t86 = -t136 * t233 + t232;
t85 = t136 * t234 + t231;
t73 = -qJD(4) * t106 - qJDD(4) * t113;
t72 = -qJD(4) * t105 + qJDD(4) * t114;
t33 = -t149 * t74 - t153 * t48;
t16 = pkin(5) * t48 - pkin(9) * t47 + t256;
t9 = t38 * qJD(5) + t160 * t49 - t164 * t50;
t4 = t163 * t5;
t1 = [qJDD(1), t191, g(1) * t166 + g(2) * t162 (t191 + (t156 ^ 2 + t158 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t176 * t157, -t176 * t155, t111 * t226 + t186 - t193, t119 * t137 - g(1) * (-pkin(1) * t162 - pkin(2) * t145 + qJ(3) * t147) - g(2) * (pkin(1) * t166 + pkin(2) * t147 + qJ(3) * t145) + t186 * t133 - t185 * qJD(3), -t104 * t105 + t114 * t70, t103 * t105 - t104 * t106 - t113 * t70 - t114 * t71, t72, t73, 0, t169 * qJD(4) + t199 * qJDD(4) + t101 * t106 + t99 * t113 + t120 * t71 + t192 * t146, -t175 * qJD(4) - t244 * qJDD(4) - t101 * t105 + t99 * t114 + t120 * t70 - t192 * t144, t182 * t47 + t34 * t75, -t182 * t48 - t34 * t74 - t35 * t75 - t47 * t62, t187, t33, 0, t136 * t192 - t37 * t149 - t9 * t153 + t62 * t256 + t80 * t35 + t69 * t48 + t53 * t74, -t135 * t192 - t38 * t149 - t8 * t153 + t182 * t256 + t80 * t34 + t69 * t47 + t53 * t75, -t58 * t212 + (t14 * t75 + t47 * t58) * t163 (-t163 * t56 - t242) * t47 + (-t12 - t15 * t163 + (t159 * t56 - t163 * t58) * qJD(6)) * t75, t247 - t259, t172 + t190, t276 * t48 + t32 * t74, -g(1) * t86 - g(2) * t88 + t37 * t15 + t4 * t74 + t6 * t48 + t9 * t56 + (t16 * t276 + t251 + (-t21 * t74 - t276 * t38 + t252) * qJD(6)) * t163 + t257 * t159, -g(1) * t85 - g(2) * t87 + t37 * t14 - t7 * t48 + t9 * t58 + (-(-qJD(6) * t38 + t16) * t276 - t251 + t206 * t74 - qJD(6) * t252) * t159 + t257 * t163; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, t155 * t84 + t157 * t83 - g(3), 0, 0, 0, 0, 0, t73, -t72, 0, 0, 0, 0, 0, t33, -t187, 0, 0, 0, 0, 0, t172 - t190, t247 + t259; 0, 0, 0, 0, -t216, t217, -t226 * qJD(1) ^ 2, t185 * qJD(1) - t264, 0, 0, 0, 0, 0, 0.2e1 * t104 * qJD(4) + t188 (-t103 - t209) * qJD(4) + t211, 0, 0, 0, 0, 0, t35 + t236, t34 - t243, 0, 0, 0, 0, 0, t181 - t250, -t163 * t276 ^ 2 - t249 - t29; 0, 0, 0, 0, 0, 0, 0, 0, t104 * t103, -t103 ^ 2 + t104 ^ 2 (t103 - t209) * qJD(4) + t211, -t188, qJDD(4), -g(3) * t146 - t101 * t104 + t193 * t144 + t201, g(3) * t144 + t101 * t103 + t193 * t146 - t184, t277, t272, t284, t274, t149, t26 * t153 + (-t104 * t62 + t149 * t164 - t153 * t223) * pkin(4) + t261, t27 * t153 + (-t104 * t182 - t149 * t160 - t153 * t222) * pkin(4) + t271, t282, -t242 * t276 + t283, t281, t181 + t250, -t266, t139 * t15 + t194 * t56 + (-t279 - t268) * t163 + t173 * t159 + t263, t139 * t14 + t194 * t58 + t173 * t163 + (-t285 + t268) * t159 + t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, t272, t284, t274, t149, t25 * t153 + t261, t24 * t153 + t271, t282, -t197 * t58 + t283, t281, -t197 * t276 + t250 + t30, -t266, -pkin(5) * t15 - t25 * t56 + t178 * t159 + (-t279 - t267) * t163 + t263, -pkin(5) * t14 - t25 * t58 + t178 * t163 + (-t285 + t267) * t159 + t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t56, -t56 ^ 2 + t58 ^ 2, t276 * t56 + t14, t276 * t58 - t15, t32, -g(1) * t87 + g(2) * t85 + t159 * t180 - t20 * t58 - t21 * t220 + t276 * t7 + t4, g(1) * t88 - g(2) * t86 + t159 * t206 + t163 * t180 + t20 * t56 + t276 * t6;];
tau_reg  = t1;
