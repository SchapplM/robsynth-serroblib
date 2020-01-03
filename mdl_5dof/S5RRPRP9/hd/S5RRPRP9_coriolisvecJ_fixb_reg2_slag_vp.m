% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:28
% EndTime: 2019-12-31 20:07:37
% DurationCPUTime: 2.83s
% Computational Cost: add. (3905->347), mult. (10012->477), div. (0->0), fcn. (6895->6), ass. (0->189)
t161 = cos(qJ(2));
t218 = t161 * qJD(1);
t144 = -qJD(4) + t218;
t157 = sin(pkin(8));
t160 = sin(qJ(2));
t224 = qJD(1) * t160;
t208 = t157 * t224;
t158 = cos(pkin(8));
t219 = t158 * qJD(2);
t112 = -t208 + t219;
t207 = t158 * t224;
t220 = t157 * qJD(2);
t113 = t207 + t220;
t159 = sin(qJ(4));
t255 = cos(qJ(4));
t65 = -t112 * t255 + t113 * t159;
t244 = t65 * t144;
t213 = t255 * t158;
t195 = t161 * t213;
t184 = qJD(1) * t195;
t217 = qJD(1) * qJD(2);
t205 = t161 * t217;
t193 = t157 * t205;
t206 = qJD(4) * t255;
t35 = t159 * (qJD(4) * t113 + t193) - qJD(2) * t184 - t112 * t206;
t266 = t35 - t244;
t265 = t65 ^ 2;
t212 = t157 * t218;
t221 = qJD(4) * t159;
t262 = -t157 * t221 + t158 * t206;
t264 = t159 * t212 - t184 + t262;
t117 = t157 * t255 + t158 * t159;
t103 = t117 * qJD(4);
t168 = t117 * t161;
t241 = -qJD(1) * t168 + t103;
t175 = t112 * t159 + t113 * t255;
t256 = t175 ^ 2;
t263 = -0.2e1 * t217;
t232 = t158 * t161;
t181 = pkin(3) * t160 - pkin(7) * t232;
t187 = pkin(2) * t160 - qJ(3) * t161;
t119 = t187 * qJD(1);
t81 = pkin(6) * t208 + t119 * t158;
t58 = qJD(1) * t181 + t81;
t104 = t157 * t119;
t233 = t158 * t160;
t234 = t157 * t161;
t171 = -pkin(6) * t233 - pkin(7) * t234;
t70 = qJD(1) * t171 + t104;
t252 = pkin(7) + qJ(3);
t128 = t252 * t157;
t129 = t252 * t158;
t79 = -t128 * t159 + t129 * t255;
t249 = qJD(3) * t117 + qJD(4) * t79 - t159 * t70 + t255 * t58;
t243 = t175 * t144;
t173 = -t159 * t157 + t213;
t155 = t160 ^ 2;
t156 = t161 ^ 2;
t261 = qJD(1) * (t155 - 0.2e1 * t156);
t260 = qJD(4) * t175;
t223 = qJD(2) * t160;
t167 = qJD(2) * t168;
t36 = qJD(1) * t167 + t260;
t54 = t160 * t262 + t167;
t95 = t117 * t160;
t259 = -t144 * t54 - t161 * t36 + (qJD(1) * t95 + t65) * t223;
t124 = -pkin(2) * t161 - qJ(3) * t160 - pkin(1);
t111 = t158 * t124;
t71 = -pkin(7) * t233 + t111 + (-pkin(6) * t157 - pkin(3)) * t161;
t235 = t157 * t160;
t142 = pkin(6) * t232;
t86 = t124 * t157 + t142;
t77 = -pkin(7) * t235 + t86;
t246 = t159 * t71 + t255 * t77;
t170 = t181 * qJD(2);
t215 = pkin(6) * t223;
t99 = qJD(2) * t187 - t160 * qJD(3);
t75 = t157 * t215 + t158 * t99;
t51 = t170 + t75;
t91 = t157 * t99;
t59 = qJD(2) * t171 + t91;
t11 = -qJD(4) * t246 - t159 * t59 + t255 * t51;
t258 = t241 * t144 - (-qJD(2) * t173 - t65) * t224;
t257 = t117 * t36 + t173 * t35 + t175 * t241 + t264 * t65;
t149 = pkin(6) * t224;
t245 = qJD(2) * pkin(2);
t201 = qJD(3) - t245;
t123 = t149 + t201;
t80 = -pkin(3) * t112 + t123;
t20 = pkin(4) * t65 - qJ(5) * t175 + t80;
t254 = t20 * t175;
t253 = t175 * t65;
t26 = t159 * t58 + t255 * t70;
t22 = qJ(5) * t224 + t26;
t174 = -t128 * t255 - t129 * t159;
t49 = qJD(3) * t173 + qJD(4) * t174;
t251 = t22 - t49;
t250 = pkin(4) * t224 + t249;
t248 = -t26 + t49;
t150 = pkin(6) * t218;
t107 = pkin(3) * t212 + t150;
t247 = -pkin(4) * t241 + qJ(5) * t264 + t117 * qJD(5) + t107;
t122 = (qJD(3) - t149) * qJD(2);
t89 = t99 * qJD(1);
t57 = t122 * t158 + t157 * t89;
t106 = t124 * qJD(1);
t131 = qJD(2) * qJ(3) + t150;
t73 = t106 * t157 + t131 * t158;
t239 = qJD(2) * t174;
t238 = qJD(2) * t79;
t237 = t112 * t158;
t163 = qJD(1) ^ 2;
t236 = t156 * t163;
t230 = t161 * t163;
t162 = qJD(2) ^ 2;
t229 = t162 * t160;
t228 = t162 * t161;
t72 = t106 * t158 - t131 * t157;
t42 = -pkin(3) * t218 - pkin(7) * t113 + t72;
t47 = pkin(7) * t112 + t73;
t17 = -t159 * t47 + t255 * t42;
t227 = qJD(5) - t17;
t143 = pkin(6) * t205;
t98 = pkin(3) * t193 + t143;
t222 = qJD(2) * t161;
t151 = pkin(6) * t222;
t211 = t161 * t220;
t108 = pkin(3) * t211 + t151;
t120 = pkin(3) * t235 + pkin(6) * t160;
t225 = t155 - t156;
t216 = pkin(6) * t234;
t214 = -t256 + t265;
t147 = -pkin(3) * t158 - pkin(2);
t209 = t144 * t224;
t148 = t160 * t217;
t56 = -t157 * t122 + t158 * t89;
t39 = qJD(1) * t170 + t56;
t46 = -pkin(7) * t193 + t57;
t204 = t159 * t46 + t206 * t47 + t221 * t42 - t255 * t39;
t202 = pkin(1) * t263;
t199 = -t112 + t219;
t198 = -t113 + t220;
t197 = pkin(4) * t148;
t196 = t174 * t35 - t36 * t79 - t49 * t65;
t192 = qJ(5) * t148;
t191 = t161 * t148;
t190 = t36 * t95 + t54 * t65;
t189 = -t256 - t265;
t188 = -t123 + t201;
t183 = qJD(1) * t199;
t182 = qJD(1) * t198;
t18 = t159 * t42 + t255 * t47;
t180 = -t144 * t18 - t204;
t30 = -t159 * t77 + t255 * t71;
t177 = -t159 * t39 - t206 * t42 + t221 * t47 - t255 * t46;
t10 = t159 * t51 + t206 * t71 - t221 * t77 + t255 * t59;
t176 = -t173 * t36 + t241 * t65;
t172 = t161 * t182;
t2 = -t197 + t204;
t53 = -qJD(2) * t195 + t103 * t160 + t159 * t211;
t96 = t173 * t160;
t166 = -t175 * t54 + t35 * t95 - t36 * t96 + t53 * t65;
t6 = pkin(4) * t36 + qJ(5) * t35 - qJD(5) * t175 + t98;
t164 = t36 - t243;
t154 = t158 ^ 2;
t153 = t157 ^ 2;
t140 = t160 * t230;
t132 = -0.2e1 * t191;
t90 = (-t144 - t218) * t223;
t85 = t111 - t216;
t82 = -pkin(6) * t207 + t104;
t76 = -t158 * t215 + t91;
t63 = -pkin(4) * t173 - qJ(5) * t117 + t147;
t43 = pkin(4) * t95 - qJ(5) * t96 + t120;
t29 = pkin(4) * t175 + qJ(5) * t65;
t28 = pkin(4) * t161 - t30;
t27 = -qJ(5) * t161 + t246;
t19 = -t35 - t244;
t16 = -t264 * t144 + (qJD(2) * t117 - t175) * t224;
t15 = -qJ(5) * t144 + t18;
t14 = pkin(4) * t144 + t227;
t13 = pkin(4) * t54 + qJ(5) * t53 - qJD(5) * t96 + t108;
t12 = -t175 * t53 - t35 * t96;
t9 = -pkin(4) * t223 - t11;
t8 = qJ(5) * t223 - qJD(5) * t161 + t10;
t7 = t53 * t144 + t35 * t161 + (qJD(1) * t96 + t175) * t223;
t5 = -t117 * t35 + t175 * t264;
t1 = -qJD(5) * t144 - t177 + t192;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t191, t225 * t263, t228, t132, -t229, 0, -pkin(6) * t228 + t160 * t202, pkin(6) * t229 + t161 * t202, 0, 0, (t113 * t158 + t154 * t224) * t222, (t237 + (-t113 - 0.2e1 * t207) * t157) * t222, (t113 * t160 + t158 * t261) * qJD(2), (-t112 * t157 + t153 * t224) * t222, (t112 * t160 - t157 * t261) * qJD(2), t132, (-qJD(1) * t75 - t56) * t161 + ((-pkin(6) * t112 + t123 * t157) * t161 + (t72 + (t85 + 0.2e1 * t216) * qJD(1)) * t160) * qJD(2), (qJD(1) * t76 + t57) * t161 + ((pkin(6) * t113 + t123 * t158) * t161 + (-t73 + (-t86 + 0.2e1 * t142) * qJD(1)) * t160) * qJD(2), t76 * t112 - t75 * t113 + (-t157 * t57 - t158 * t56) * t160 + (-t157 * t73 - t158 * t72 + (-t157 * t86 - t158 * t85) * qJD(1)) * t222, t56 * t85 + t57 * t86 + t72 * t75 + t73 * t76 + (t123 + t149) * t151, t12, t166, t7, t190, -t259, t90, t108 * t65 - t11 * t144 + t120 * t36 + t204 * t161 + t80 * t54 + t98 * t95 + (qJD(1) * t30 + t17) * t223, t10 * t144 + t108 * t175 - t120 * t35 - t177 * t161 - t80 * t53 + t98 * t96 + (-qJD(1) * t246 - t18) * t223, -t10 * t65 - t11 * t175 + t17 * t53 + t177 * t95 - t18 * t54 + t204 * t96 - t246 * t36 + t30 * t35, t10 * t18 + t108 * t80 + t11 * t17 + t120 * t98 - t177 * t246 - t204 * t30, t12, t7, -t166, t90, t259, t190, t13 * t65 + t9 * t144 + t2 * t161 + t20 * t54 + t43 * t36 + t6 * t95 + (-qJD(1) * t28 - t14) * t223, -t1 * t95 - t14 * t53 - t15 * t54 + t175 * t9 + t2 * t96 - t27 * t36 - t28 * t35 - t65 * t8, -t1 * t161 - t13 * t175 - t8 * t144 + t20 * t53 + t43 * t35 - t6 * t96 + (qJD(1) * t27 + t15) * t223, t1 * t27 + t13 * t20 + t14 * t9 + t15 * t8 + t2 * t28 + t43 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, t225 * t163, 0, t140, 0, 0, t163 * pkin(1) * t160, pkin(1) * t230, 0, 0, t158 * t172, (-t237 + t113 * t157 + (-t153 + t154) * qJD(2)) * t218, t158 * t236 + t160 * t182, -t199 * t212, -t157 * t236 + t160 * t183, t140, ((-qJ(3) * t220 - t72) * t160 + (-pkin(6) * t199 + t157 * t188 + t81) * t161) * qJD(1), ((-qJ(3) * t219 + t73) * t160 + (pkin(6) * t198 + t158 * t188 - t82) * t161) * qJD(1), -t82 * t112 + t81 * t113 + (qJD(3) * t112 + t218 * t72 + t57) * t158 + (qJD(3) * t113 + t218 * t73 - t56) * t157, -t72 * t81 - t73 * t82 + (-t157 * t72 + t158 * t73) * qJD(3) + (-t56 * t157 + t57 * t158) * qJ(3) + (-t123 - t245) * t150, t5, -t257, t16, t176, t258, t209, -t107 * t65 - t98 * t173 + t147 * t36 + t241 * t80 + t249 * t144 + (-t17 + t239) * t224, -t107 * t175 + t98 * t117 - t147 * t35 + t264 * t80 + t248 * t144 + (t18 - t238) * t224, t117 * t204 - t17 * t264 - t173 * t177 + t175 * t249 - t18 * t241 + t26 * t65 + t196, -t80 * t107 + t98 * t147 - t17 * t249 - t174 * t204 - t177 * t79 + t18 * t248, t5, t16, t257, t209, -t258, t176, -t6 * t173 + t63 * t36 - t247 * t65 + t241 * t20 + t250 * t144 + (t14 + t239) * t224, t1 * t173 + t2 * t117 + t14 * t264 - t15 * t241 + t175 * t250 + t22 * t65 + t196, -t6 * t117 + t63 * t35 + t247 * t175 - t264 * t20 + t251 * t144 + (-t15 + t238) * t224, t1 * t79 + t14 * t250 - t15 * t251 - t174 * t2 - t20 * t247 + t6 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t161 * t183, -t112 ^ 2 - t113 ^ 2, -t112 * t73 + t113 * t72 + t143, 0, 0, 0, 0, 0, 0, t164, -t266, t189, t17 * t175 + t18 * t65 + t98, 0, 0, 0, 0, 0, 0, t164, t189, t266, -t14 * t175 + t15 * t65 + t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, -t214, t19, -t253, -t117 * t205 - t243 - t260, t148, -t175 * t80 + t180, -t144 * t17 + t65 * t80 + t177, 0, 0, t253, t19, t214, t148, t36 + t243, -t253, -t29 * t65 + t180 + 0.2e1 * t197 - t254, pkin(4) * t35 - t36 * qJ(5) + (t15 - t18) * t175 + (t14 - t227) * t65, 0.2e1 * t192 - t20 * t65 + t29 * t175 + (-0.2e1 * qJD(5) + t17) * t144 - t177, -t2 * pkin(4) + t1 * qJ(5) - t14 * t18 + t15 * t227 - t20 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148 + t253, t19, -t144 ^ 2 - t256, t144 * t15 + t2 + t254;];
tauc_reg = t3;
