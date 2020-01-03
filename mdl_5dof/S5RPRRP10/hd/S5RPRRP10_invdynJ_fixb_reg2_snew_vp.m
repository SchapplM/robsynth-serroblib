% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP10
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:04
% EndTime: 2019-12-31 18:52:12
% DurationCPUTime: 2.86s
% Computational Cost: add. (9022->317), mult. (21642->426), div. (0->0), fcn. (15807->8), ass. (0->208)
t176 = sin(pkin(8));
t173 = t176 ^ 2;
t177 = cos(pkin(8));
t174 = t177 ^ 2;
t221 = t173 + t174;
t181 = sin(qJ(1));
t248 = cos(qJ(1));
t193 = t248 * g(1) + t181 * g(2);
t276 = (2 * qJD(2) * qJD(1)) - t193;
t183 = cos(qJ(3));
t180 = sin(qJ(3));
t229 = t177 * t180;
t195 = t176 * t183 + t229;
t163 = t195 * qJD(1);
t179 = sin(qJ(4));
t182 = cos(qJ(4));
t147 = -t182 * qJD(3) + t179 * t163;
t149 = t179 * qJD(3) + t182 * t163;
t117 = t149 * t147;
t213 = t177 * qJDD(1);
t214 = t176 * qJDD(1);
t197 = t180 * t214 - t183 * t213;
t217 = t163 * qJD(3);
t135 = -t197 - t217;
t127 = qJDD(4) - t135;
t264 = -t117 + t127;
t275 = pkin(4) * t264;
t160 = t195 * qJDD(1);
t161 = (t176 * t180 - t177 * t183) * qJD(1);
t218 = t161 * qJD(3);
t137 = t160 - t218;
t108 = -t147 * qJD(4) + t179 * qJDD(3) + t182 * t137;
t156 = qJD(4) + t161;
t124 = t156 * t147;
t88 = t108 + t124;
t274 = qJ(5) * t88;
t273 = t179 * t264;
t139 = t163 * t161;
t263 = qJDD(3) - t139;
t272 = t180 * t263;
t271 = t182 * t264;
t270 = t183 * t263;
t184 = qJD(1) ^ 2;
t269 = -(t184 * pkin(1)) + qJDD(1) * qJ(2) + t276;
t146 = t149 ^ 2;
t155 = t156 ^ 2;
t112 = -t146 - t155;
t145 = t147 ^ 2;
t201 = -t182 * qJDD(3) + t179 * t137;
t107 = -t149 * qJD(4) - t201;
t119 = t156 * pkin(4) - t149 * qJ(5);
t128 = t161 * pkin(3) - t163 * pkin(7);
t207 = pkin(2) * t177 + pkin(1);
t266 = qJ(2) + pkin(6);
t188 = (-qJDD(1) * t266 + t207 * t184 - t276) * t176;
t246 = t177 * g(3);
t187 = t188 - t246;
t199 = -t176 * g(3) + t269 * t177;
t123 = -t174 * t184 * pkin(2) + pkin(6) * t213 + t199;
t226 = t183 * t123;
t259 = qJD(3) ^ 2;
t70 = -t259 * pkin(3) + qJDD(3) * pkin(7) - t161 * t128 + t180 * t187 + t226;
t203 = t181 * g(1) - t248 * g(2);
t198 = -qJDD(2) + t203;
t130 = t207 * qJDD(1) + (t221 * pkin(6) + qJ(2)) * t184 + t198;
t73 = (-t137 + t218) * pkin(7) + (-t135 + t217) * pkin(3) - t130;
t43 = t179 * t73 + t182 * t70;
t192 = t107 * qJ(5) - 0.2e1 * qJD(5) * t147 - t156 * t119 + t43;
t268 = -t192 + (t112 + t145) * pkin(4);
t265 = t108 - t124;
t223 = t184 * qJ(2);
t232 = qJDD(1) * pkin(1);
t157 = t198 + t223 + t232;
t262 = t221 * t223 - t157 - t232;
t85 = (qJD(4) - t156) * t149 + t201;
t158 = t161 ^ 2;
t159 = t163 ^ 2;
t109 = -t155 - t145;
t62 = t179 * t109 + t271;
t258 = pkin(3) * t62;
t96 = t117 + t127;
t238 = t179 * t96;
t67 = t182 * t112 - t238;
t257 = pkin(3) * t67;
t219 = qJD(5) * t149;
t142 = -0.2e1 * t219;
t42 = t179 * t70 - t182 * t73;
t191 = -t274 - t42 + t275;
t24 = t142 + t191;
t256 = pkin(4) * t24;
t255 = pkin(4) * t88;
t103 = -t145 - t146;
t56 = t179 * t88 - t182 * t85;
t38 = -t183 * t103 + t180 * t56;
t254 = pkin(6) * t38;
t63 = t182 * t109 - t273;
t84 = (qJD(4) + t156) * t149 + t201;
t45 = t180 * t63 - t183 * t84;
t253 = pkin(6) * t45;
t236 = t182 * t96;
t68 = -t179 * t112 - t236;
t49 = t180 * t68 - t183 * t265;
t252 = pkin(6) * t49;
t54 = -t179 * t85 - t182 * t88;
t251 = pkin(7) * t54;
t250 = pkin(7) * t62;
t249 = pkin(7) * t67;
t247 = pkin(3) * t180;
t39 = t180 * t103 + t183 * t56;
t245 = qJ(2) * (-t176 * t38 + t177 * t39) - pkin(1) * t54;
t46 = t180 * t84 + t183 * t63;
t244 = qJ(2) * (-t176 * t45 + t177 * t46) - pkin(1) * t62;
t50 = t180 * t265 + t183 * t68;
t243 = qJ(2) * (-t176 * t49 + t177 * t50) - pkin(1) * t67;
t242 = -pkin(3) * t84 + pkin(7) * t63;
t241 = -pkin(3) * t265 + pkin(7) * t68;
t91 = t180 * t123 - t183 * t187;
t92 = -g(3) * t229 + t180 * t188 + t226;
t58 = t180 * t92 - t183 * t91;
t240 = t176 * t58;
t69 = -qJDD(3) * pkin(3) - t259 * pkin(7) + t163 * t128 + t91;
t239 = t179 * t69;
t237 = t182 * t69;
t235 = -pkin(3) * t103 + pkin(7) * t56;
t234 = qJ(5) * t179;
t233 = qJ(5) * t182;
t231 = t156 * t179;
t230 = t156 * t182;
t228 = t180 * t130;
t132 = qJDD(3) + t139;
t227 = t180 * t132;
t225 = t183 * t130;
t224 = t183 * t132;
t209 = t180 * t117;
t208 = t183 * t117;
t206 = -pkin(3) * t183 - pkin(2);
t205 = -pkin(2) * t62 + pkin(6) * t46;
t204 = -pkin(2) * t67 + pkin(6) * t50;
t20 = t179 * t42 + t182 * t43;
t59 = t180 * t91 + t183 * t92;
t200 = t176 * (t269 * t176 + t246) + t177 * t199;
t19 = t179 * t43 - t182 * t42;
t190 = t191 + t275;
t35 = -t107 * pkin(4) - t145 * qJ(5) + t149 * t119 + qJDD(5) + t69;
t169 = t174 * qJDD(1);
t168 = t173 * qJDD(1);
t164 = t221 * t184;
t152 = -t159 - t259;
t151 = -t159 + t259;
t150 = t158 - t259;
t143 = 0.2e1 * t219;
t136 = t160 - 0.2e1 * t218;
t134 = t197 + 0.2e1 * t217;
t129 = -t259 - t158;
t121 = -t146 + t155;
t120 = t145 - t155;
t116 = -t158 - t159;
t114 = t146 - t145;
t111 = -t180 * t152 - t224;
t110 = t183 * t152 - t227;
t102 = t180 * t160 - t183 * t197;
t101 = -t183 * t160 - t180 * t197;
t99 = t183 * t129 - t272;
t98 = t180 * t129 + t270;
t94 = (-t147 * t182 + t149 * t179) * t156;
t93 = (-t147 * t179 - t149 * t182) * t156;
t81 = t182 * t108 - t149 * t231;
t80 = t179 * t108 + t149 * t230;
t79 = -t179 * t107 + t147 * t230;
t78 = t182 * t107 + t147 * t231;
t77 = t182 * t120 - t238;
t76 = -t179 * t121 + t271;
t75 = t179 * t120 + t236;
t74 = t182 * t121 + t273;
t57 = -pkin(4) * t265 - qJ(5) * t96;
t55 = -t179 * t265 - t182 * t84;
t53 = -t179 * t84 + t182 * t265;
t47 = t176 * (t180 * t127 + t183 * t94) + t177 * (-t183 * t127 + t180 * t94);
t40 = t237 - t249;
t37 = pkin(6) * t39;
t36 = t239 - t250;
t34 = -pkin(3) * t54 + t255;
t33 = t176 * (t183 * t81 + t209) + t177 * (t180 * t81 - t208);
t32 = t176 * (t183 * t79 - t209) + t177 * (t180 * t79 + t208);
t31 = -qJ(5) * t112 + t35;
t30 = t43 - t257;
t29 = t42 - t258;
t28 = -t145 * pkin(4) + t192;
t27 = -pkin(4) * t84 + qJ(5) * t109 - t35;
t26 = t176 * (-t180 * t85 + t183 * t77) + t177 * (t180 * t77 + t183 * t85);
t25 = t176 * (t180 * t88 + t183 * t76) + t177 * (t180 * t76 - t183 * t88);
t21 = t176 * (t180 * t114 + t183 * t55) + t177 * (-t183 * t114 + t180 * t55);
t17 = t143 - t191 + t274;
t16 = -qJ(5) * t85 + (-t103 - t145) * pkin(4) + t192;
t15 = -t257 - t268;
t14 = -t179 * t57 + t182 * t31 - t249;
t13 = -t179 * t27 - t233 * t264 - t250;
t10 = t143 - t190 - t258;
t9 = -pkin(4) * t35 + qJ(5) * t28;
t8 = -t19 - t251;
t7 = -t179 * t24 + t182 * t28;
t6 = t179 * t28 + t182 * t24;
t5 = t180 * t35 + t183 * t7;
t4 = t180 * t7 - t183 * t35;
t3 = -pkin(3) * t6 - t256;
t2 = -t179 * t16 + t182 * t17 - t251;
t1 = -pkin(7) * t6 - t179 * t9 - t24 * t233;
t11 = [0, 0, 0, 0, 0, qJDD(1), t203, t193, 0, 0, t168, 0.2e1 * t176 * t213, 0, t169, 0, 0, -t262 * t177, t262 * t176, pkin(1) * t164 + qJ(2) * (t169 + t168) + t200, pkin(1) * t157 + qJ(2) * t200, t176 * (t183 * t137 - t180 * t217) + t177 * (t180 * t137 + t183 * t217), t176 * (-t183 * t134 - t180 * t136) + t177 * (-t180 * t134 + t183 * t136), t176 * (-t180 * t151 + t270) + t177 * (t183 * t151 + t272), t176 * (-t180 * t135 + t183 * t218) + t177 * (t183 * t135 + t180 * t218), t176 * (t183 * t150 - t227) + t177 * (t180 * t150 + t224), (t176 * (-t161 * t183 + t163 * t180) + t177 * (-t161 * t180 - t163 * t183)) * qJD(3), t176 * (-pkin(6) * t98 - t228) + t177 * (-pkin(2) * t134 + pkin(6) * t99 + t225) - pkin(1) * t134 + qJ(2) * (-t176 * t98 + t177 * t99), t176 * (-pkin(6) * t110 - t225) + t177 * (-pkin(2) * t136 + pkin(6) * t111 - t228) - pkin(1) * t136 + qJ(2) * (-t176 * t110 + t177 * t111), t176 * (-pkin(6) * t101 - t58) + t177 * (-pkin(2) * t116 + pkin(6) * t102 + t59) - pkin(1) * t116 + qJ(2) * (-t176 * t101 + t177 * t102), -pkin(6) * t240 + t177 * (pkin(2) * t130 + pkin(6) * t59) + pkin(1) * t130 + qJ(2) * (t177 * t59 - t240), t33, t21, t25, t32, t26, t47, t176 * (-t180 * t29 + t183 * t36 - t253) + t177 * (t180 * t36 + t183 * t29 + t205) + t244, t176 * (-t180 * t30 + t183 * t40 - t252) + t177 * (t180 * t40 + t183 * t30 + t204) + t243, t176 * (t183 * t8 + t247 * t54 - t254) + t177 * (t180 * t8 + t206 * t54 + t37) + t245, (t176 * (-pkin(7) * t183 + t247) + t177 * (-pkin(7) * t180 + t206) - pkin(1)) * t19 + t266 * (-t176 * (t180 * t20 - t183 * t69) + t177 * (t180 * t69 + t183 * t20)), t33, t21, t25, t32, t26, t47, t176 * (-t180 * t10 + t183 * t13 - t253) + t177 * (t183 * t10 + t180 * t13 + t205) + t244, t176 * (t183 * t14 - t180 * t15 - t252) + t177 * (t180 * t14 + t183 * t15 + t204) + t243, t176 * (-t180 * t34 + t183 * t2 - t254) + t177 * (-pkin(2) * t54 + t180 * t2 + t183 * t34 + t37) + t245, t176 * (-pkin(6) * t4 + t183 * t1 - t180 * t3) + t177 * (-pkin(2) * t6 + pkin(6) * t5 + t180 * t1 + t183 * t3) - pkin(1) * t6 + qJ(2) * (-t176 * t4 + t177 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213, t214, -t164, -t157, 0, 0, 0, 0, 0, 0, t134, t136, t116, -t130, 0, 0, 0, 0, 0, 0, t62, t67, t54, t19, 0, 0, 0, 0, 0, 0, t62, t67, t54, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t159 - t158, t160, -t139, -t197, qJDD(3), -t91, -t92, 0, 0, t80, t53, t74, t78, t75, t93, -t237 + t242, t239 + t241, t20 + t235, -pkin(3) * t69 + pkin(7) * t20, t80, t53, t74, t78, t75, t93, t182 * t27 - t234 * t264 + t242, t179 * t31 + t182 * t57 + t241, t182 * t16 + t179 * t17 + t235, -pkin(3) * t35 + pkin(7) * t7 + t182 * t9 - t234 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, t114, t88, -t117, -t85, t127, -t42, -t43, 0, 0, t117, t114, t88, -t117, -t85, t127, t142 + t190, t268, -t255, t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t265, t103, t35;];
tauJ_reg = t11;
