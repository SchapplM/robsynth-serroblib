% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:34
% EndTime: 2019-03-09 05:13:41
% DurationCPUTime: 2.73s
% Computational Cost: add. (5124->311), mult. (13633->398), div. (0->0), fcn. (10873->8), ass. (0->175)
t128 = qJD(3) + qJD(4);
t132 = sin(qJ(6));
t134 = cos(qJ(6));
t130 = sin(pkin(10));
t131 = cos(pkin(10));
t219 = sin(qJ(3));
t221 = cos(qJ(3));
t153 = t130 * t219 - t131 * t221;
t108 = t153 * qJD(1);
t152 = -t130 * t221 - t131 * t219;
t109 = t152 * qJD(1);
t133 = sin(qJ(4));
t220 = cos(qJ(4));
t87 = t108 * t220 - t109 * t133;
t76 = t128 * t132 - t134 * t87;
t211 = t76 * t87;
t158 = -t108 * t133 - t109 * t220;
t228 = qJD(6) + t158;
t245 = t134 * t228;
t143 = t220 * t153;
t140 = qJD(3) * t143;
t181 = qJD(4) * t220;
t192 = qJD(4) * t133;
t182 = qJD(3) * t219;
t183 = qJD(3) * t221;
t195 = (t130 * t183 + t131 * t182) * qJD(1);
t50 = qJD(1) * t140 + t108 * t181 - t109 * t192 + t133 * t195;
t252 = t132 * t50 - t228 * t245;
t259 = t252 - t211;
t78 = t128 * t134 + t132 * t87;
t214 = t78 * t87;
t246 = t132 * t228;
t254 = -t134 * t50 - t228 * t246;
t258 = t254 + t214;
t207 = t128 * t87;
t151 = -t50 + t207;
t212 = t87 ^ 2;
t239 = t158 ^ 2;
t257 = -t212 + t239;
t256 = qJD(6) - t228;
t190 = qJD(6) * t134;
t191 = qJD(6) * t132;
t142 = t108 * qJD(3);
t51 = -t108 * t192 - t109 * t181 - t133 * t142 + t195 * t220;
t28 = -t128 * t191 + t132 * t51 + t190 * t87;
t27 = t28 * t134;
t255 = -t246 * t78 + t27;
t250 = t228 * t76;
t49 = t134 * t51;
t29 = qJD(6) * t78 - t49;
t253 = -t245 * t78 + (-t28 + t250) * t132 - t134 * t29;
t223 = pkin(5) * t87;
t210 = t158 * t87;
t209 = pkin(7) + qJ(2);
t115 = t209 * t130;
t112 = qJD(1) * t115;
t116 = t209 * t131;
t113 = qJD(1) * t116;
t75 = -pkin(8) * t108 - t112 * t219 + t113 * t221;
t72 = t220 * t75;
t229 = -t112 * t221 - t113 * t219;
t74 = pkin(8) * t109 + t229;
t73 = qJD(3) * pkin(3) + t74;
t39 = t133 * t73 + t72;
t37 = -qJ(5) * t128 - t39;
t20 = -t37 - t223;
t251 = t20 * t228;
t249 = t228 * t87;
t202 = t87 * qJ(5);
t59 = -pkin(8) * t195 - qJD(2) * t108 + qJD(3) * t229;
t184 = qJD(2) * t219;
t185 = qJD(2) * t221;
t230 = -t130 * t185 - t131 * t184;
t60 = pkin(8) * t142 + qJD(1) * t230 + t112 * t182 - t113 * t183;
t13 = t133 * t59 + t181 * t75 + t192 * t73 - t220 * t60;
t123 = -t131 * pkin(2) - pkin(1);
t114 = qJD(1) * t123 + qJD(2);
t94 = t108 * pkin(3) + t114;
t141 = -qJ(5) * t158 + t94;
t40 = pkin(4) * t87 + t141;
t216 = t40 * t158;
t248 = t13 + t216;
t71 = t133 * t75;
t38 = -t220 * t73 + t71;
t198 = qJD(5) + t38;
t179 = -t133 * t60 - t181 * t73 + t192 * t75 - t220 * t59;
t12 = -qJD(5) * t128 + t179;
t3 = -pkin(5) * t51 - t12;
t237 = pkin(5) * t158;
t199 = t237 + t198;
t224 = pkin(4) + pkin(9);
t18 = -t128 * t224 + t199;
t23 = t224 * t87 + t141;
t8 = t132 * t18 + t134 * t23;
t244 = t134 * t3 - t8 * t87;
t243 = -t40 * t87 - t12;
t242 = t87 * t94 + t179;
t165 = t132 * t23 - t134 * t18;
t241 = t132 * t3 - t165 * t87 + t190 * t20;
t238 = pkin(4) * t158;
t236 = t158 * t20;
t235 = t158 * t94;
t42 = t220 * t74 - t71;
t200 = pkin(3) * t181 + qJD(5) - t42;
t41 = t133 * t74 + t72;
t169 = pkin(3) * t192 - t41;
t203 = t158 * t128;
t233 = -t51 + t203;
t232 = t158 * t224;
t124 = -pkin(3) * t220 - pkin(4);
t121 = -pkin(9) + t124;
t227 = (t169 + t223) * t228 - t121 * t50;
t225 = qJD(3) ^ 2;
t92 = -t133 * t152 + t143;
t218 = t20 * t92;
t148 = t133 * t153;
t93 = -t152 * t220 - t148;
t98 = pkin(3) * t153 + t123;
t138 = -t93 * qJ(5) + t98;
t31 = t224 * t92 + t138;
t217 = t31 * t50;
t215 = t50 * t92;
t139 = -t115 * t183 - t116 * t182 - t130 * t184 + t131 * t185;
t146 = qJD(3) * t152;
t63 = pkin(8) * t146 + t139;
t147 = qJD(3) * t153;
t64 = pkin(8) * t147 + t115 * t182 - t116 * t183 + t230;
t83 = pkin(8) * t152 - t115 * t221 - t116 * t219;
t154 = t115 * t219 - t116 * t221;
t84 = -pkin(8) * t153 - t154;
t15 = -t133 * t64 - t181 * t83 + t192 * t84 - t220 * t63;
t205 = t15 * t128;
t159 = t133 * t83 + t220 * t84;
t16 = qJD(4) * t159 + t133 * t63 - t220 * t64;
t204 = t16 * t128;
t201 = t200 + t237;
t194 = t130 ^ 2 + t131 ^ 2;
t193 = qJD(3) * t109;
t189 = qJD(1) * qJD(2);
t187 = t92 * t190;
t180 = t195 * pkin(3);
t174 = t194 * qJD(1) ^ 2;
t45 = t133 * t84 - t220 * t83;
t58 = -qJD(4) * t148 - t133 * t147 - t146 * t220 - t152 * t181;
t168 = t228 * t58 - t215;
t167 = -t109 * pkin(3) + t202;
t166 = -t224 * t50 + (t39 - t223) * t228;
t164 = (-qJD(6) * t121 + t167 + t232) * t228;
t163 = (qJD(6) * t224 + t202 + t232) * t228;
t162 = 0.2e1 * t194 * t189;
t160 = t128 * t39 - t13;
t32 = pkin(5) * t93 + t45;
t157 = t20 * t58 + t3 * t92 + t32 * t50;
t150 = t50 + t207;
t149 = t50 * qJ(5) - qJD(5) * t158 + t180;
t145 = t152 * qJD(2);
t144 = pkin(3) * t146;
t14 = t51 * pkin(4) + t149;
t57 = qJD(4) * t143 - t133 * t146 - t152 * t192 + t140;
t17 = t58 * pkin(4) + t57 * qJ(5) - t93 * qJD(5) - t144;
t137 = t51 + t203;
t122 = pkin(3) * t133 + qJ(5);
t52 = t202 + t238;
t44 = t92 * pkin(4) + t138;
t43 = t167 + t238;
t36 = -pkin(4) * t128 + t198;
t35 = t50 * t93;
t33 = -pkin(5) * t92 + t159;
t11 = t58 * pkin(9) + t17;
t10 = -t57 * pkin(5) + t16;
t9 = -pkin(5) * t58 - t15;
t6 = t224 * t51 + t149;
t5 = -pkin(5) * t50 + t13;
t4 = t134 * t5;
t1 = [0, 0, 0, 0, 0, t162, qJ(2) * t162 (t108 * t152 + t109 * t153) * qJD(3), qJD(1) * qJD(3) * t153 ^ 2 + t108 * t147 - t109 * t146 + t152 * t195, -t153 * t225, t152 * t225, 0, t123 * t195 + (qJD(3) * t154 - t114 * t152 + t145) * qJD(3), -qJD(3) * t139 - t114 * t147 - t123 * t142, -t158 * t57 - t35, -t158 * t58 - t51 * t93 + t57 * t87 + t215, -t57 * t128, -t58 * t128, 0, -t144 * t87 + t180 * t92 + t98 * t51 + t94 * t58 - t204, -t144 * t158 + t180 * t93 - t98 * t50 - t94 * t57 + t205, t12 * t92 + t13 * t93 + t15 * t87 + t158 * t16 - t159 * t51 - t36 * t57 + t37 * t58 - t45 * t50, -t14 * t92 - t17 * t87 - t40 * t58 - t44 * t51 + t204, -t14 * t93 - t158 * t17 + t40 * t57 + t44 * t50 - t205, -t12 * t159 + t13 * t45 + t14 * t44 + t15 * t37 + t16 * t36 + t17 * t40, t78 * t187 + (t28 * t92 + t58 * t78) * t132 (-t132 * t76 + t134 * t78) * t58 + (-t132 * t29 + t27 + (-t132 * t78 - t134 * t76) * qJD(6)) * t92, t132 * t168 + t187 * t228 + t28 * t93 - t57 * t78, -t191 * t228 * t92 + t134 * t168 - t29 * t93 + t57 * t76, -t228 * t57 - t35, t33 * t29 + t4 * t93 + t165 * t57 + t9 * t76 + (-t11 * t228 - t6 * t93 + t217) * t132 + (t10 * t228 - t157) * t134 + ((-t132 * t32 - t134 * t31) * t228 - t8 * t93 + t132 * t218) * qJD(6), t33 * t28 + t8 * t57 + t9 * t78 + (-(qJD(6) * t32 + t11) * t228 + t217 - (qJD(6) * t18 + t6) * t93 + qJD(6) * t218) * t134 + (-(-qJD(6) * t31 + t10) * t228 - (-qJD(6) * t23 + t5) * t93 + t157) * t132; 0, 0, 0, 0, 0, -t174, -qJ(2) * t174, 0, 0, 0, 0, 0, -t193 + t195, -0.2e1 * t142, 0, 0, 0, 0, 0, t137, -t150, -t212 - t239, -t137, t150, -t158 * t36 - t37 * t87 + t14, 0, 0, 0, 0, 0, t252 + t211, t214 - t254; 0, 0, 0, 0, 0, 0, 0, -t109 * t108, -t108 ^ 2 + t109 ^ 2, 0, -t193 - t195, 0, qJD(1) * t145 + t114 * t109, t114 * t108 + t153 * t189, t210, t257, t151, t233, 0, t128 * t41 - t235 + (t109 * t87 - t128 * t192) * pkin(3) - t13, t42 * t128 + (t109 * t158 - t128 * t181) * pkin(3) + t242, -t122 * t51 - t124 * t50 + (t169 - t37) * t158 + (-t200 + t36) * t87, t128 * t169 + t43 * t87 + t248, t128 * t200 + t158 * t43 + t243, -t12 * t122 + t124 * t13 + t169 * t36 - t200 * t37 - t40 * t43, t255, t253, t258, t259, t249, t122 * t29 + t201 * t76 + t132 * t164 + (t227 + t236) * t134 + t241, t122 * t28 + t201 * t78 + t134 * t164 + (-t227 - t251) * t132 + t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, t257, t151, t233, 0, t160 - t235, -t128 * t38 + t242, pkin(4) * t50 - qJ(5) * t51 + (-t37 - t39) * t158 + (t36 - t198) * t87, t52 * t87 - t160 + t216, t128 * t198 + t158 * t52 + t243, -pkin(4) * t13 - qJ(5) * t12 - t198 * t37 - t36 * t39 - t40 * t52, t255, t253, t258, t259, t249, qJ(5) * t29 + t199 * t76 + t132 * t163 + (-t166 + t236) * t134 + t241, qJ(5) * t28 + t199 * t78 + t134 * t163 + (t166 - t251) * t132 + t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, -t210, -t128 ^ 2 - t239, t128 * t37 + t248, 0, 0, 0, 0, 0, -t128 * t76 + t254, -t128 * t78 + t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 * t76, -t76 ^ 2 + t78 ^ 2, t28 + t250, -t256 * t78 + t49, -t50, -t132 * t6 - t20 * t78 - t256 * t8 + t4, -t132 * t5 - t134 * t6 + t165 * t256 + t20 * t76;];
tauc_reg  = t1;
