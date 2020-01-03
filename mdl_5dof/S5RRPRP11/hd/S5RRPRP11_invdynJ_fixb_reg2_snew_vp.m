% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRP11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:04
% EndTime: 2019-12-31 20:14:11
% DurationCPUTime: 2.29s
% Computational Cost: add. (4122->268), mult. (8584->287), div. (0->0), fcn. (4822->6), ass. (0->180)
t150 = cos(qJ(4));
t147 = sin(qJ(4));
t151 = cos(qJ(2));
t195 = qJD(1) * qJD(2);
t135 = t151 * t195;
t148 = sin(qJ(2));
t137 = t148 * qJDD(1);
t114 = t137 + t135;
t104 = qJDD(4) + t114;
t198 = qJD(1) * t151;
t108 = t147 * qJD(2) + t150 * t198;
t110 = t150 * qJD(2) - t147 * t198;
t211 = t110 * t108;
t60 = -t211 - t104;
t226 = t147 * t60;
t103 = t110 ^ 2;
t197 = t148 * qJD(1);
t132 = qJD(4) + t197;
t238 = t132 ^ 2;
t247 = -t103 - t238;
t29 = t150 * t247 + t226;
t281 = pkin(3) * t29;
t280 = t148 * t29;
t236 = -pkin(2) - pkin(7);
t279 = t236 * t29;
t207 = t148 * qJ(3);
t165 = t151 * t236 - pkin(1) - t207;
t217 = t150 * t60;
t278 = t165 * (-t147 * t247 + t217);
t186 = t148 * t195;
t192 = t151 * qJDD(1);
t115 = -t186 + t192;
t181 = t147 * qJDD(2) + t150 * t115;
t167 = t110 * qJD(4) + t181;
t91 = t132 * t110;
t40 = t167 - t91;
t239 = t108 ^ 2;
t85 = t239 - t238;
t277 = -t148 * t40 + t151 * (-t147 * t85 + t217);
t212 = t108 * t132;
t170 = t150 * qJDD(2) - t147 * t115;
t65 = -t108 * qJD(4) + t170;
t250 = t65 - t212;
t276 = qJ(5) * t250;
t273 = -t150 * t85 - t226;
t219 = t150 * t250;
t39 = t167 + t91;
t73 = t103 - t239;
t272 = t151 * (-t147 * t39 + t219) - t148 * t73;
t54 = -t239 - t103;
t271 = pkin(3) * t54;
t244 = -t211 + t104;
t216 = t150 * t244;
t243 = -t238 - t239;
t253 = t147 * t243 + t216;
t270 = pkin(3) * t253;
t269 = qJ(3) * t54;
t265 = t151 * t54;
t264 = t148 * t253;
t263 = t236 * t253;
t87 = -t103 + t238;
t262 = -t147 * t87 + t216;
t249 = t65 + t212;
t51 = t147 * t244;
t260 = t148 * t249 + t151 * (-t150 * t87 - t51);
t259 = 2 * qJD(3);
t143 = t148 ^ 2;
t154 = qJD(1) ^ 2;
t138 = t143 * t154;
t153 = qJD(2) ^ 2;
t126 = -t138 - t153;
t202 = t148 * t154;
t188 = t151 * t202;
t121 = qJDD(2) - t188;
t201 = t151 * t121;
t258 = pkin(6) * (t148 * t126 + t201);
t252 = t147 * t250 + t150 * t39;
t251 = t165 * (t150 * t243 - t51);
t235 = pkin(3) + pkin(6);
t174 = t114 + t135;
t248 = t174 * qJ(3);
t231 = t151 * pkin(2);
t175 = -t207 - t231;
t111 = t175 * qJD(1);
t140 = t148 * g(3);
t149 = sin(qJ(1));
t152 = cos(qJ(1));
t177 = t152 * g(1) + t149 * g(2);
t213 = qJDD(1) * pkin(6);
t96 = -t154 * pkin(1) - t177 + t213;
t163 = (qJD(1) * t111 + t96) * t151 - t153 * pkin(2) - t140;
t161 = qJD(2) * t259 + t163;
t144 = t151 ^ 2;
t208 = t144 * t154;
t246 = t201 + t148 * (-t153 + t208);
t242 = pkin(4) * t167 - t276;
t193 = qJDD(2) * qJ(3);
t241 = t115 * pkin(3) - pkin(7) * t208 + t193;
t116 = -0.2e1 * t186 + t192;
t128 = -t153 - t208;
t120 = qJDD(2) + t188;
t205 = t148 * t120;
t240 = pkin(6) * (-t151 * t128 + t205) - pkin(1) * t116;
t237 = 2 * qJD(5);
t222 = t148 * t96;
t166 = -qJDD(2) * pkin(2) - t153 * qJ(3) + t111 * t197 + qJDD(3) + t222;
t155 = t114 * pkin(3) - qJDD(2) * pkin(7) + (-pkin(3) * t195 - pkin(7) * t202 + g(3)) * t151 + t166;
t123 = pkin(3) * t197 - qJD(2) * pkin(7);
t131 = pkin(2) * t186;
t187 = qJD(3) * t197;
t134 = -0.2e1 * t187;
t185 = t149 * g(1) - t152 * g(2);
t169 = -qJDD(1) * pkin(1) - t185;
t21 = -t123 * t197 + t131 + t134 + (-pkin(3) * t144 - pkin(6)) * t154 + t236 * t115 - t248 + t169;
t15 = t147 * t155 + t150 * t21;
t70 = t108 * pkin(4) - t110 * qJ(5);
t180 = t104 * qJ(5) - t108 * t70 + t132 * t237 + t15;
t10 = -pkin(4) * t238 + t180;
t14 = t147 * t21 - t150 * t155;
t12 = -t104 * pkin(4) - qJ(5) * t238 + t110 * t70 + qJDD(5) + t14;
t234 = -pkin(4) * t12 + qJ(5) * t10;
t232 = pkin(4) * t147;
t230 = t151 * g(3);
t229 = -pkin(4) * t249 - qJ(5) * t40;
t28 = (t259 + t123) * qJD(2) + t163 + t241;
t228 = t147 * t28;
t218 = t150 * t249;
t42 = (-qJD(4) + t132) * t110 - t181;
t17 = t147 * t42 - t218;
t224 = t148 * t17;
t221 = t150 * t28;
t214 = qJ(5) * t147;
t210 = t132 * t147;
t209 = t132 * t150;
t206 = t148 * t116;
t118 = t138 + t208;
t200 = pkin(1) * t118 + (t143 + t144) * t213;
t196 = qJD(4) + t132;
t191 = t108 * t210;
t190 = t148 * t211;
t189 = t108 * t209;
t81 = t222 + t230;
t82 = t151 * t96 - t140;
t184 = t148 * t81 + t151 * t82;
t183 = -qJ(5) * t150 + qJ(3);
t84 = t110 * t209;
t179 = t151 * (-t147 * t65 - t84) + t190;
t83 = t110 * t210;
t178 = t83 - t189;
t3 = -t150 * t14 + t147 * t15;
t173 = t147 * t14 + t150 * t15;
t172 = t151 * (-t138 + t153) + t205;
t168 = t147 * t167 + t189;
t164 = t148 * t104 + t151 * (t84 + t191);
t95 = t154 * pkin(6) - t169;
t162 = -pkin(4) * t247 - qJ(5) * t60 + t10;
t49 = t166 + t230;
t160 = pkin(4) * t244 + qJ(5) * t243 - t12;
t159 = t115 * pkin(2) - t131 + t95;
t158 = -t190 + t151 * (t150 * t167 - t191);
t48 = t161 + t193;
t157 = t159 + 0.2e1 * t187;
t156 = -qJD(2) * t123 + t110 * t237 - t161 - t241 - t242;
t119 = t138 - t208;
t113 = t137 + 0.2e1 * t135;
t80 = t174 * t148;
t79 = (t115 - t186) * t151;
t68 = t151 * t113 + t206;
t45 = -t196 * t108 + t170;
t41 = t196 * t110 + t181;
t36 = t147 * t249;
t35 = t150 * t65 - t83;
t16 = -t147 * t40 - t218;
t13 = (pkin(4) * t132 - (2 * qJD(5))) * t110 + t28 + t242;
t8 = t156 + (-t41 - t91) * pkin(4);
t7 = -pkin(4) * t91 + t156 + t276;
t6 = -qJ(5) * t54 + t12;
t5 = (-t238 - t54) * pkin(4) + t180;
t1 = t147 * t10 - t150 * t12;
t2 = [0, 0, 0, 0, 0, qJDD(1), t185, t177, 0, 0, t80, t68, t172, t79, t246, 0, t151 * t95 - t240, -pkin(1) * t113 - t148 * t95 - t258, t184 + t200, pkin(1) * t95 + pkin(6) * t184, 0, -t172, -t246, t80, t68, t79, t148 * (qJ(3) * t118 + t166) + (pkin(2) * t118 + t140 + t48) * t151 + t200, t151 * (-pkin(2) * t116 + t134 - t159) + (-t151 * t174 - t206) * qJ(3) + t240, t148 * t157 + t258 + (pkin(1) + t231) * t113 + (t113 + t174) * t207, pkin(6) * (t148 * t49 + t151 * t48) + (pkin(1) - t175) * (t157 + t248), t179, -t272, t260, t158, t277, t164, t148 * (-t14 + t270) + t151 * (pkin(3) * t39 + t221) + pkin(6) * (t151 * t39 + t264) + t251, t148 * (-t15 + t281) + t151 * (pkin(3) * t45 - t228) + pkin(6) * (t151 * t45 + t280) + t278, pkin(3) * t224 + t151 * (-t173 + t271) + pkin(6) * (t224 + t265) + t165 * (t150 * t42 + t36), t165 * t173 + t235 * (t148 * t3 + t151 * t28), t179, t260, t272, t164, -t277, t158, t148 * (t160 + t270) + t151 * (-t150 * t8 + (pkin(3) + t214) * t41) + pkin(6) * (t151 * t41 + t264) + t251, t148 * (pkin(3) * t16 + t229) + t151 * (-t147 * t6 - t150 * t5 + t271) + pkin(6) * (t148 * t16 + t265) + t165 * (-t150 * t40 + t36), t148 * (t162 - t281) + t151 * (-pkin(3) * t250 - pkin(4) * t219 - t147 * t7) + pkin(6) * (-t151 * t250 - t280) - t278, (t235 * t1 + t234) * t148 + t165 * (t150 * t10 + t147 * t12) + (pkin(4) * t150 + t214 + t235) * t151 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, t119, t137, t188, t192, qJDD(2), -t81, -t82, 0, 0, qJDD(2), -t137, -t192, -t188, t119, t188, (-pkin(2) * t148 + qJ(3) * t151) * qJDD(1), -pkin(2) * t120 - qJ(3) * t128 + t49, -pkin(2) * t126 + (qJDD(2) + t121) * qJ(3) + t161, -pkin(2) * t49 + qJ(3) * t48, t35, -t252, t262, t168, -t273, t178, qJ(3) * t39 + t228 + t263, qJ(3) * t45 + t221 + t279, t236 * t17 + t269 - t3, qJ(3) * t28 + t236 * t3, t35, t262, t252, t178, t273, t168, -t147 * t8 + t183 * t41 + t263, -t147 * t5 + t150 * t6 + t236 * t16 + t269, t150 * t7 + (-qJ(3) - t232) * t250 - t279, t236 * t1 + (t183 + t232) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t120, t126, t49, 0, 0, 0, 0, 0, 0, t253, t29, t17, t3, 0, 0, 0, 0, 0, 0, t253, t16, -t29, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t73, t249, -t211, -t40, t104, -t14, -t15, 0, 0, t211, t249, -t73, t104, t40, -t211, t160, t229, t162, t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, t249, t247, t12;];
tauJ_reg = t2;
