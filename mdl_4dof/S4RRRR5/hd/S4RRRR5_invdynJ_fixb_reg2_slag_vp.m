% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:16
% EndTime: 2019-12-31 17:28:24
% DurationCPUTime: 3.72s
% Computational Cost: add. (3502->390), mult. (8133->544), div. (0->0), fcn. (5349->10), ass. (0->194)
t150 = cos(qJ(2));
t204 = t150 * qJD(1);
t123 = -qJD(3) + t204;
t146 = sin(qJ(2));
t147 = sin(qJ(1));
t151 = cos(qJ(1));
t179 = g(1) * t151 + g(2) * t147;
t248 = g(3) * t150;
t162 = t146 * t179 - t248;
t201 = t146 * qJDD(1);
t130 = pkin(5) * t201;
t203 = qJD(1) * qJD(2);
t189 = t150 * t203;
t230 = qJDD(2) * pkin(2);
t74 = pkin(5) * t189 + t130 - t230;
t159 = -t74 + t162;
t270 = -qJD(3) * pkin(6) * t123 - t159;
t149 = cos(qJ(3));
t223 = t149 * t150;
t171 = pkin(3) * t146 - pkin(7) * t223;
t257 = pkin(7) + pkin(6);
t196 = qJD(3) * t257;
t145 = sin(qJ(3));
t214 = qJD(1) * t146;
t195 = t145 * t214;
t180 = pkin(2) * t146 - pkin(6) * t150;
t97 = t180 * qJD(1);
t55 = pkin(5) * t195 + t149 * t97;
t269 = -qJD(1) * t171 - t149 * t196 - t55;
t226 = t146 * t149;
t227 = t145 * t150;
t81 = t145 * t97;
t268 = t81 + (-pkin(5) * t226 - pkin(7) * t227) * qJD(1) + t145 * t196;
t144 = sin(qJ(4));
t148 = cos(qJ(4));
t205 = t149 * qJD(2);
t92 = t195 - t205;
t191 = t149 * t214;
t206 = t145 * qJD(2);
t94 = t191 + t206;
t175 = t144 * t92 - t148 * t94;
t39 = t144 * t94 + t148 * t92;
t247 = t39 * t175;
t267 = qJD(2) * qJD(3) + t189 + t201;
t266 = t175 ^ 2 - t39 ^ 2;
t115 = -qJD(4) + t123;
t188 = qJD(3) * t214;
t184 = t267 * t145 + t149 * t188;
t166 = t149 * qJDD(2) - t184;
t208 = qJD(4) * t148;
t209 = qJD(4) * t144;
t37 = (-qJDD(2) + t188) * t145 - t267 * t149;
t8 = -t144 * t166 + t148 * t37 + t92 * t208 + t94 * t209;
t265 = -t115 * t39 - t8;
t132 = pkin(5) * t204;
t107 = qJD(2) * pkin(6) + t132;
t181 = t150 * pkin(2) + t146 * pkin(6);
t102 = -pkin(1) - t181;
t86 = t102 * qJD(1);
t47 = -t107 * t145 + t149 * t86;
t31 = -t94 * pkin(7) + t47;
t27 = -t123 * pkin(3) + t31;
t48 = t149 * t107 + t145 * t86;
t32 = -t92 * pkin(7) + t48;
t100 = t180 * qJD(2);
t51 = qJD(1) * t100 + qJDD(1) * t102;
t44 = t149 * t51;
t134 = t150 * qJDD(1);
t187 = t146 * t203;
t260 = -t187 + t134;
t73 = t260 * pkin(5) + qJDD(2) * pkin(6);
t15 = -qJD(3) * t48 - t145 * t73 + t44;
t89 = qJDD(3) - t260;
t6 = t89 * pkin(3) + t37 * pkin(7) + t15;
t210 = qJD(3) * t149;
t211 = qJD(3) * t145;
t14 = -t107 * t211 + t145 * t51 + t149 * t73 + t86 * t210;
t7 = pkin(7) * t166 + t14;
t1 = (qJD(4) * t27 + t7) * t148 + t144 * t6 - t32 * t209;
t143 = qJ(3) + qJ(4);
t137 = cos(t143);
t249 = g(3) * t146;
t106 = -qJD(2) * pkin(2) + pkin(5) * t214;
t57 = pkin(3) * t92 + t106;
t136 = sin(t143);
t221 = t151 * t136;
t224 = t147 * t150;
t66 = -t137 * t224 + t221;
t220 = t151 * t137;
t68 = t147 * t136 + t150 * t220;
t264 = g(1) * t68 - g(2) * t66 + t137 * t249 + t57 * t39 - t1;
t240 = t148 * t32;
t11 = t144 * t27 + t240;
t2 = -qJD(4) * t11 - t144 * t7 + t148 * t6;
t65 = t136 * t224 + t220;
t67 = t147 * t137 - t150 * t221;
t263 = -g(1) * t67 + g(2) * t65 + t136 * t249 + t57 * t175 + t2;
t158 = qJD(4) * t175 + t144 * t37 + t148 * t166;
t262 = t115 * t175 + t158;
t261 = t47 * t123 + t14;
t129 = t149 * pkin(3) + pkin(2);
t173 = t150 * t129 + t146 * t257;
t218 = t151 * t149;
t77 = t145 * t224 + t218;
t219 = t151 * t145;
t79 = t147 * t149 - t150 * t219;
t259 = -g(1) * t79 + g(2) * t77;
t200 = qJD(3) + qJD(4);
t258 = -0.2e1 * pkin(1);
t254 = pkin(3) * t145;
t253 = g(1) * t147;
t250 = g(2) * t151;
t246 = t94 * t92;
t108 = t257 * t145;
t109 = t257 * t149;
t53 = -t148 * t108 - t144 * t109;
t245 = qJD(4) * t53 + t269 * t144 - t268 * t148;
t54 = -t144 * t108 + t148 * t109;
t244 = -qJD(4) * t54 + t268 * t144 + t269 * t148;
t229 = t144 * t145;
t95 = -t148 * t149 + t229;
t243 = -t148 * t210 - t149 * t208 + t200 * t229 - t95 * t204;
t96 = t144 * t149 + t148 * t145;
t50 = t200 * t96;
t242 = -t96 * t204 + t50;
t241 = t144 * t32;
t239 = t37 * t145;
t237 = t48 * t123;
t236 = t92 * t123;
t235 = t92 * t145;
t234 = t94 * t123;
t233 = t94 * t149;
t232 = t146 * pkin(5) * t206 + t149 * t100;
t124 = pkin(5) * t223;
t60 = t145 * t102 + t124;
t231 = pkin(5) * qJDD(1);
t228 = t145 * t146;
t217 = t151 * pkin(1) + t147 * pkin(5);
t141 = t146 ^ 2;
t142 = t150 ^ 2;
t216 = t141 - t142;
t215 = t141 + t142;
t213 = qJD(2) * t146;
t212 = qJD(2) * t150;
t207 = t106 * qJD(3);
t154 = qJD(1) ^ 2;
t197 = t146 * t154 * t150;
t194 = t150 * t205;
t193 = t146 * t210;
t192 = t150 * t206;
t183 = t150 * t187;
t182 = pkin(3) * t211 - t204 * t254 - t132;
t178 = pkin(5) * t92 + t106 * t145;
t177 = pkin(5) * t94 + t106 * t149;
t176 = -pkin(6) * t89 + t207;
t91 = t149 * t102;
t46 = -pkin(7) * t226 + t91 + (-pkin(5) * t145 - pkin(3)) * t150;
t52 = -pkin(7) * t228 + t60;
t19 = -t144 * t52 + t148 * t46;
t20 = t144 * t46 + t148 * t52;
t174 = -t145 * t48 - t149 * t47;
t170 = -pkin(5) * qJDD(2) + t203 * t258;
t169 = -t123 * t210 + t145 * t89;
t168 = t123 * t211 + t149 * t89;
t167 = pkin(1) * t154 + t179;
t153 = qJD(2) ^ 2;
t165 = pkin(5) * t153 + qJDD(1) * t258 + t250;
t164 = t192 + t193;
t163 = t166 * t149;
t160 = -t150 * t179 - t249;
t29 = t145 * t100 + t102 * t210 + (-t146 * t205 - t150 * t211) * pkin(5);
t139 = t151 * pkin(5);
t125 = t146 * t253;
t101 = (pkin(5) + t254) * t146;
t85 = qJDD(4) + t89;
t80 = t147 * t145 + t150 * t218;
t78 = -t147 * t223 + t219;
t71 = t95 * t146;
t70 = t96 * t146;
t59 = -pkin(5) * t227 + t91;
t58 = pkin(3) * t164 + pkin(5) * t212;
t56 = -pkin(5) * t191 + t81;
t30 = -qJD(3) * t60 + t232;
t26 = -pkin(3) * t166 + t74;
t23 = (t200 * t226 + t192) * t148 + (-t200 * t228 + t194) * t144;
t22 = t144 * t192 + t146 * t50 - t148 * t194;
t21 = -pkin(7) * t164 + t29;
t18 = t171 * qJD(2) + (-t124 + (pkin(7) * t146 - t102) * t145) * qJD(3) + t232;
t13 = t148 * t31 - t241;
t12 = -t144 * t31 - t240;
t10 = t148 * t27 - t241;
t4 = -qJD(4) * t20 - t144 * t21 + t148 * t18;
t3 = qJD(4) * t19 + t144 * t18 + t148 * t21;
t5 = [0, 0, 0, 0, 0, qJDD(1), -t250 + t253, t179, 0, 0, t141 * qJDD(1) + 0.2e1 * t183, 0.2e1 * t134 * t146 - 0.2e1 * t203 * t216, qJDD(2) * t146 + t153 * t150, t142 * qJDD(1) - 0.2e1 * t183, qJDD(2) * t150 - t153 * t146, 0, t170 * t146 + (-t165 + t253) * t150, t146 * t165 + t150 * t170 - t125, 0.2e1 * t215 * t231 - t179, -g(1) * (-t147 * pkin(1) + t139) - g(2) * t217 + (pkin(5) ^ 2 * t215 + pkin(1) ^ 2) * qJDD(1), t94 * t194 + (-t37 * t149 - t211 * t94) * t146, (-t94 * t145 - t149 * t92) * t212 + (t163 + t239 + (-t233 + t235) * qJD(3)) * t146, (-t123 * t205 + t37) * t150 + (qJD(2) * t94 + t168) * t146, t92 * t193 + (-t146 * t166 + t212 * t92) * t145, (t123 * t206 - t166) * t150 + (-t92 * qJD(2) - t169) * t146, -t123 * t213 - t89 * t150, -g(1) * t78 - g(2) * t80 - t30 * t123 + t59 * t89 + (qJD(2) * t178 - t15) * t150 + (-pkin(5) * t166 + t47 * qJD(2) + t74 * t145 + t149 * t207) * t146, -g(1) * t77 - g(2) * t79 + t29 * t123 - t60 * t89 + (qJD(2) * t177 + t14) * t150 + (-pkin(5) * t37 - qJD(2) * t48 - t145 * t207 + t74 * t149) * t146, -t29 * t92 + t60 * t166 - t30 * t94 + t59 * t37 + t125 + t174 * t212 + (-t250 - t14 * t145 - t15 * t149 + (t145 * t47 - t149 * t48) * qJD(3)) * t146, t14 * t60 + t48 * t29 + t15 * t59 + t47 * t30 - g(1) * t139 - g(2) * (t151 * t181 + t217) - t102 * t253 + (t106 * t212 + t74 * t146) * pkin(5), t175 * t22 + t71 * t8, -t158 * t71 + t175 * t23 + t22 * t39 + t70 * t8, t22 * t115 + t8 * t150 - t175 * t213 - t71 * t85, -t158 * t70 + t23 * t39, t23 * t115 - t150 * t158 - t213 * t39 - t70 * t85, -t115 * t213 - t85 * t150, -g(1) * t66 - g(2) * t68 + t10 * t213 - t101 * t158 - t4 * t115 - t2 * t150 + t19 * t85 + t57 * t23 + t26 * t70 + t58 * t39, -g(1) * t65 - g(2) * t67 + t1 * t150 - t101 * t8 - t11 * t213 + t3 * t115 - t175 * t58 - t20 * t85 - t57 * t22 - t26 * t71, -t1 * t70 + t10 * t22 - t11 * t23 - t146 * t250 + t158 * t20 + t175 * t4 + t19 * t8 + t2 * t71 - t3 * t39 + t125, t1 * t20 + t11 * t3 + t2 * t19 + t10 * t4 + t26 * t101 + t57 * t58 - g(1) * (pkin(3) * t219 + t139) - g(2) * (t173 * t151 + t217) + (-g(1) * (-pkin(1) - t173) - g(2) * t254) * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t197, t216 * t154, t201, t197, t134, qJDD(2), t146 * t167 - t130 - t248, t249 + (t167 - t231) * t150, 0, 0, -t123 * t233 - t239, (-t37 + t236) * t149 + (t166 + t234) * t145, (t123 * t223 - t146 * t94) * qJD(1) + t169, -t123 * t235 + t163, (-t123 * t227 + t92 * t146) * qJD(1) + t168, t123 * t214, -pkin(2) * t184 + t55 * t123 + t176 * t145 + (-t47 * t146 - t150 * t178) * qJD(1) + (t230 - t270) * t149, pkin(2) * t37 - t56 * t123 + t176 * t149 + (t146 * t48 - t150 * t177) * qJD(1) + t270 * t145, t55 * t94 + t56 * t92 + ((t94 * qJD(3) + t166) * pkin(6) + t261) * t149 + (-t15 + t237 + (t92 * qJD(3) - t37) * pkin(6)) * t145 + t160, -t106 * t132 - t47 * t55 - t48 * t56 + t159 * pkin(2) + (qJD(3) * t174 + t14 * t149 - t15 * t145 + t160) * pkin(6), t175 * t243 - t8 * t96, t158 * t96 + t175 * t242 + t243 * t39 + t8 * t95, t243 * t115 + t175 * t214 + t96 * t85, -t158 * t95 + t242 * t39, t242 * t115 + t39 * t214 - t95 * t85, t115 * t214, -t10 * t214 - t244 * t115 + t129 * t158 + t162 * t137 + t182 * t39 + t242 * t57 + t26 * t95 + t53 * t85, t11 * t214 + t245 * t115 + t129 * t8 - t136 * t162 - t175 * t182 - t243 * t57 + t26 * t96 - t54 * t85, -t1 * t95 + t243 * t10 - t242 * t11 + t158 * t54 + t175 * t244 - t2 * t96 - t245 * t39 + t53 * t8 + t160, -g(3) * t173 + t1 * t54 + t244 * t10 + t245 * t11 - t26 * t129 + t182 * t57 + t2 * t53 + t179 * (t129 * t146 - t150 * t257); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, -t92 ^ 2 + t94 ^ 2, -t37 - t236, -t246, t166 - t234, t89, -t107 * t210 - t106 * t94 - t237 + t44 + (-qJD(3) * t86 + t249 - t73) * t145 + t259, g(1) * t80 - g(2) * t78 + g(3) * t226 + t106 * t92 - t261, 0, 0, -t247, t266, t265, t247, t262, t85, t12 * t115 + (t115 * t209 + t148 * t85 - t39 * t94) * pkin(3) + t263, -t13 * t115 + (t115 * t208 - t144 * t85 + t175 * t94) * pkin(3) + t264, -t10 * t39 - t11 * t175 - t12 * t175 + t13 * t39 + (t144 * t158 + t148 * t8 + (-t144 * t175 - t148 * t39) * qJD(4)) * pkin(3), -t10 * t12 - t11 * t13 + (t1 * t144 + t2 * t148 - t57 * t94 + g(3) * t228 + (-t10 * t144 + t11 * t148) * qJD(4) + t259) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t247, t266, t265, t247, t262, t85, -t11 * t115 + t263, -t10 * t115 + t264, 0, 0;];
tau_reg = t5;
