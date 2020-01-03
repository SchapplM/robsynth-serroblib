% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:10
% EndTime: 2019-12-31 18:52:16
% DurationCPUTime: 3.48s
% Computational Cost: add. (4703->402), mult. (11327->485), div. (0->0), fcn. (8403->10), ass. (0->210)
t150 = sin(pkin(8));
t155 = sin(qJ(3));
t151 = cos(pkin(8));
t265 = cos(qJ(3));
t206 = t265 * t151;
t173 = -t155 * t150 + t206;
t280 = t173 * qJD(1);
t281 = t280 * qJD(3);
t232 = qJDD(1) * pkin(1);
t140 = qJDD(2) - t232;
t156 = sin(qJ(1));
t158 = cos(qJ(1));
t279 = -g(1) * t156 + g(2) * t158;
t177 = -t140 - t279;
t117 = t265 * t150 + t155 * t151;
t169 = t117 * qJDD(1);
t162 = t169 + t281;
t278 = qJD(3) * qJD(4) + t162;
t111 = t117 * qJD(1);
t154 = sin(qJ(4));
t157 = cos(qJ(4));
t214 = qJD(4) * t154;
t171 = t154 * qJDD(3) - t111 * t214 + t278 * t157;
t250 = qJ(5) * t171;
t114 = t117 * qJD(3);
t210 = t150 * qJDD(1);
t184 = -qJDD(1) * t206 + t155 * t210;
t79 = qJD(1) * t114 + t184;
t73 = qJDD(4) + t79;
t267 = pkin(4) * t73;
t137 = t151 * pkin(2) + pkin(1);
t123 = -t137 * qJD(1) + qJD(2);
t57 = -pkin(3) * t280 - t111 * pkin(7) + t123;
t254 = pkin(6) + qJ(2);
t124 = t254 * t150;
t118 = qJD(1) * t124;
t125 = t254 * t151;
t119 = qJD(1) * t125;
t81 = -t155 * t118 + t265 * t119;
t75 = qJD(3) * pkin(7) + t81;
t32 = t154 * t57 + t157 * t75;
t203 = qJD(3) * t265;
t212 = qJD(1) * qJD(2);
t270 = t254 * qJDD(1) + t212;
t95 = t270 * t150;
t96 = t270 * t151;
t207 = t118 * t203 + t155 * t95 - t265 * t96;
t215 = qJD(3) * t155;
t42 = -t119 * t215 - t207;
t37 = qJDD(3) * pkin(7) + t42;
t209 = t151 * qJDD(1);
t41 = -pkin(2) * t209 + t79 * pkin(3) - pkin(7) * t162 + t140;
t7 = -qJD(4) * t32 - t154 * t37 + t157 * t41;
t216 = qJD(3) * t154;
t89 = t111 * t157 + t216;
t1 = -qJD(5) * t89 - t250 + t267 + t7;
t100 = qJD(4) - t280;
t87 = -t157 * qJD(3) + t111 * t154;
t25 = -qJ(5) * t87 + t32;
t248 = t100 * t25;
t277 = t1 + t248;
t246 = t100 * t32;
t276 = t7 + t246;
t149 = pkin(8) + qJ(3);
t141 = sin(t149);
t189 = g(1) * t158 + g(2) * t156;
t175 = t189 * t141;
t275 = -t265 * t124 - t155 * t125;
t272 = qJ(2) * qJDD(1);
t142 = cos(t149);
t219 = t157 * t158;
t224 = t154 * t156;
t101 = t142 * t224 + t219;
t220 = t156 * t157;
t223 = t154 * t158;
t103 = -t142 * t223 + t220;
t256 = g(3) * t154;
t271 = -g(1) * t103 + g(2) * t101 + t141 * t256;
t258 = g(3) * t141;
t167 = -t189 * t142 - t258;
t269 = t89 ^ 2;
t268 = t111 ^ 2;
t266 = t87 * pkin(4);
t128 = t158 * t137;
t260 = g(2) * t128;
t257 = g(3) * t142;
t255 = t89 * t87;
t253 = qJ(5) + pkin(7);
t31 = -t154 * t75 + t157 * t57;
t24 = -qJ(5) * t89 + t31;
t23 = pkin(4) * t100 + t24;
t252 = -t24 + t23;
t213 = qJD(4) * t157;
t187 = -t157 * qJDD(3) + t111 * t213;
t51 = t154 * t169 + (qJD(4) + t280) * t216 + t187;
t251 = -t154 * t51 - t87 * t213;
t76 = pkin(3) * t111 - pkin(7) * t280;
t108 = t155 * t119;
t80 = -t265 * t118 - t108;
t45 = t154 * t76 + t157 * t80;
t78 = -pkin(3) * t173 - pkin(7) * t117 - t137;
t85 = -t155 * t124 + t265 * t125;
t82 = t157 * t85;
t49 = t154 * t78 + t82;
t249 = qJ(5) * t51;
t247 = t100 * t31;
t245 = t100 * t87;
t244 = t111 * t87;
t243 = t154 * t73;
t242 = t154 * t87;
t241 = t154 * t89;
t240 = t157 * t87;
t239 = t171 * t154;
t238 = t51 * t157;
t237 = t89 * t100;
t236 = t89 * t111;
t198 = qJD(4) * t253;
t230 = t280 * t154;
t235 = qJ(5) * t230 + t157 * qJD(5) - t154 * t198 - t45;
t44 = -t154 * t80 + t157 * t76;
t234 = -pkin(4) * t111 - t154 * qJD(5) - t44 + (qJ(5) * t280 - t198) * t157;
t233 = qJD(4) * t89;
t231 = t100 * t111;
t229 = t111 * t280;
t113 = t150 * t215 - t151 * t203;
t228 = t113 * t154;
t227 = t113 * t157;
t226 = t117 * t154;
t225 = t117 * t157;
t218 = (g(1) * t219 + g(2) * t220) * t141;
t147 = t150 ^ 2;
t148 = t151 ^ 2;
t217 = t147 + t148;
t58 = t173 * qJD(2) + t275 * qJD(3);
t77 = pkin(3) * t114 + pkin(7) * t113;
t208 = t154 * t77 + t157 * t58 + t78 * t213;
t205 = t117 * t213;
t43 = t118 * t215 - t119 * t203 - t155 * t96 - t265 * t95;
t38 = -qJDD(3) * pkin(3) - t43;
t16 = pkin(4) * t51 + qJDD(5) + t38;
t204 = -t16 - t257;
t202 = qJD(5) + t266;
t201 = pkin(4) * t154 + t254;
t199 = -t154 * t58 + t157 * t77;
t48 = -t154 * t85 + t157 * t78;
t6 = t154 * t41 + t157 * t37 + t57 * t213 - t75 * t214;
t197 = t217 * qJD(1) ^ 2;
t196 = t100 * t157;
t195 = 0.2e1 * t217;
t194 = pkin(7) * qJD(4) * t100 + t38;
t193 = t279 * t141;
t192 = pkin(3) * t142 + pkin(7) * t141;
t191 = -g(1) * t101 - g(2) * t103;
t102 = -t142 * t220 + t223;
t104 = t142 * t219 + t224;
t190 = -g(1) * t102 - g(2) * t104;
t2 = -qJD(5) * t87 - t249 + t6;
t186 = -t100 * t23 + t2;
t185 = t6 - t247;
t183 = t154 * t32 + t157 * t31;
t182 = t240 + t241;
t139 = pkin(4) * t157 + pkin(3);
t181 = t139 * t142 + t141 * t253;
t179 = qJ(5) * t113 - qJD(5) * t117;
t178 = t157 * t73 + (-t214 + t230) * t100;
t176 = t157 * t171 - t89 * t214;
t74 = -qJD(3) * pkin(3) - t80;
t174 = -pkin(7) * t73 + t100 * t74;
t172 = t205 - t228;
t122 = -t137 * qJDD(1) + qJDD(2);
t168 = t177 + t232;
t166 = -t257 + t175;
t165 = g(1) * t104 - g(2) * t102 + t157 * t258 - t6;
t164 = t195 * t212 - t189;
t163 = t7 + t271;
t59 = qJD(2) * t117 + qJD(3) * t85;
t161 = t278 * t154 + t187;
t131 = t142 * t256;
t127 = t253 * t157;
t126 = t253 * t154;
t107 = t280 ^ 2;
t86 = t87 ^ 2;
t60 = pkin(4) * t226 - t275;
t53 = pkin(4) * t230 + t81;
t52 = t202 + t74;
t46 = -t86 + t269;
t40 = pkin(4) * t172 + t59;
t39 = t100 * t114 - t173 * t73;
t33 = -qJ(5) * t226 + t49;
t29 = -pkin(4) * t173 - qJ(5) * t225 + t48;
t28 = -t161 + t237;
t27 = t171 + t245;
t22 = -t100 ^ 2 * t157 - t236 - t243;
t21 = t100 * t196 - t236 + t243;
t20 = t178 + t244;
t19 = t178 - t244;
t18 = t100 * t242 - t238;
t17 = t196 * t89 + t239;
t15 = -t49 * qJD(4) + t199;
t14 = -t85 * t214 + t208;
t13 = -t251 * t117 - t87 * t228;
t12 = t117 * t176 - t89 * t227;
t11 = -qJ(5) * t205 + (-qJD(4) * t85 + t179) * t154 + t208;
t10 = pkin(4) * t114 + t179 * t157 + (-t82 + (qJ(5) * t117 - t78) * t154) * qJD(4) + t199;
t9 = -t100 * t172 - t87 * t114 + t173 * t51 - t73 * t226;
t8 = t73 * t225 + t89 * t114 - t171 * t173 + (-t117 * t214 - t227) * t100;
t5 = t182 * t280 + t176 + t251;
t4 = -(-t240 + t241) * t280 - t176 + t251;
t3 = t182 * t113 + (-t239 - t238 + (-t157 * t89 + t242) * qJD(4)) * t117;
t26 = [0, 0, 0, 0, 0, qJDD(1), -t279, t189, 0, 0, t147 * qJDD(1), 0.2e1 * t150 * t209, 0, t148 * qJDD(1), 0, 0, t168 * t151, -t168 * t150, t195 * t272 + t164, t177 * pkin(1) + (t217 * t272 + t164) * qJ(2), -t111 * t113 + t117 * t162, -t111 * t114 - t113 * t280 - t117 * t79 + t162 * t173, -qJD(3) * t113 + qJDD(3) * t117, -t114 * t280 - t173 * t79, -qJD(3) * t114 + qJDD(3) * t173, 0, -qJD(3) * t59 + qJDD(3) * t275 + t114 * t123 - t122 * t173 - t137 * t79 - t142 * t279, -t58 * qJD(3) - t85 * qJDD(3) - t123 * t113 + t122 * t117 - t137 * t162 + t193, t59 * t111 + t80 * t113 - t81 * t114 - t43 * t117 - t162 * t275 + t42 * t173 + t280 * t58 - t85 * t79 - t189, t42 * t85 + t81 * t58 + t43 * t275 - t80 * t59 - t122 * t137 - g(1) * (-t156 * t137 + t158 * t254) - g(2) * (t156 * t254 + t128), t12, t3, t8, t13, t9, t39, -t74 * t228 + t100 * t15 + t114 * t31 - t173 * t7 + t48 * t73 - t51 * t275 + t59 * t87 + (t154 * t38 + t213 * t74) * t117 + t190, -t74 * t227 - t100 * t14 - t114 * t32 + t173 * t6 - t49 * t73 - t171 * t275 + t59 * t89 + (t157 * t38 - t214 * t74) * t117 + t191, -t14 * t87 - t15 * t89 - t48 * t171 - t49 * t51 + t183 * t113 + (-t154 * t6 - t157 * t7 + (t154 * t31 - t157 * t32) * qJD(4)) * t117 - t193, -t260 + t32 * t14 + t31 * t15 - t38 * t275 + t7 * t48 + t6 * t49 + t74 * t59 + (-g(1) * t254 - g(2) * t192) * t158 + (-g(1) * (-t137 - t192) - g(2) * t254) * t156, t12, t3, t8, t13, t9, t39, -t52 * t228 - t1 * t173 + t10 * t100 + t114 * t23 + t29 * t73 + t40 * t87 + t51 * t60 + (t154 * t16 + t213 * t52) * t117 + t190, -t52 * t227 - t100 * t11 - t114 * t25 + t173 * t2 - t33 * t73 + t40 * t89 + t171 * t60 + (t157 * t16 - t214 * t52) * t117 + t191, -t10 * t89 - t11 * t87 - t29 * t171 - t33 * t51 + (t154 * t25 + t157 * t23) * t113 + (-t1 * t157 - t154 * t2 + (t154 * t23 - t157 * t25) * qJD(4)) * t117 - t193, -t260 + t1 * t29 + t23 * t10 + t25 * t11 + t16 * t60 + t2 * t33 + t52 * t40 + (-g(1) * t201 - g(2) * t181) * t158 + (-g(1) * (-t137 - t181) - g(2) * t201) * t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t210, -t197, -qJ(2) * t197 - t177, 0, 0, 0, 0, 0, 0, 0.2e1 * t111 * qJD(3) + t184, t169 + 0.2e1 * t281, -t107 - t268, t111 * t80 - t280 * t81 + t122 + t279, 0, 0, 0, 0, 0, 0, t19, t22, t4, -t111 * t74 + t185 * t154 + t276 * t157 + t279, 0, 0, 0, 0, 0, 0, t19, t22, t4, -t111 * t52 + t186 * t154 + t277 * t157 + t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t229, -t107 + t268, t169, t229, -t184, qJDD(3), qJD(3) * t81 - t111 * t123 + t166 + t43, -t280 * t123 + (t80 + t108) * qJD(3) + t207 - t167, 0, 0, t17, t5, t21, t18, t20, -t231, -pkin(3) * t51 - t100 * t44 - t111 * t31 - t81 * t87 + (-t194 - t257) * t157 + t174 * t154 + t218, -pkin(3) * t171 + t100 * t45 + t111 * t32 - t81 * t89 + t131 + t174 * t157 + (-t175 + t194) * t154, t44 * t89 + t45 * t87 + ((-t51 + t233) * pkin(7) + t185) * t157 + ((qJD(4) * t87 + t171) * pkin(7) - t276) * t154 + t167, -t31 * t44 - t32 * t45 - t74 * t81 + (-t38 + t166) * pkin(3) + (-qJD(4) * t183 - t7 * t154 + t6 * t157 + t167) * pkin(7), t17, t5, t21, t18, t20, -t231, -t111 * t23 - t126 * t73 - t139 * t51 - t53 * t87 + t204 * t157 + t234 * t100 + (-t280 * t52 + (t52 + t266) * qJD(4)) * t154 + t218, t111 * t25 - t127 * t73 - t139 * t171 - t53 * t89 + t131 + t52 * t196 - t235 * t100 + (pkin(4) * t233 + t16 - t175) * t154, t126 * t171 - t127 * t51 - t277 * t154 + t186 * t157 - t234 * t89 - t235 * t87 + t167, t2 * t127 - t1 * t126 - t16 * t139 - g(3) * t181 + (pkin(4) * t214 - t53) * t52 + t235 * t25 + t234 * t23 + t189 * (t139 * t141 - t142 * t253); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t255, t46, t27, -t255, t28, t73, -t74 * t89 + t163 + t246, t74 * t87 + t165 + t247, 0, 0, t255, t46, t27, -t255, t28, t73, 0.2e1 * t267 - t250 + t248 + (-t202 - t52) * t89 + t163, -pkin(4) * t269 + t249 + t100 * t24 + (qJD(5) + t52) * t87 + t165, -pkin(4) * t171 - t252 * t87, t252 * t25 + (-t52 * t89 + t1 + t271) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161 + t237, t171 - t245, -t86 - t269, t23 * t89 + t25 * t87 - t175 - t204;];
tau_reg = t26;
