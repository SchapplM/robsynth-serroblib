% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:26:00
% EndTime: 2022-01-23 09:26:06
% DurationCPUTime: 2.44s
% Computational Cost: add. (4083->273), mult. (10949->398), div. (0->0), fcn. (7976->8), ass. (0->177)
t154 = cos(pkin(8));
t201 = t154 * qJD(1);
t140 = -qJD(3) + t201;
t152 = sin(pkin(8));
t151 = sin(pkin(9));
t153 = cos(pkin(9));
t156 = sin(qJ(3));
t158 = cos(qJ(3));
t127 = t151 * t158 + t153 * t156;
t167 = qJD(1) * t127;
t101 = t152 * t167;
t210 = qJD(1) * t152;
t192 = t156 * t210;
t178 = t151 * t192;
t191 = t158 * t210;
t105 = t153 * t191 - t178;
t155 = sin(qJ(5));
t157 = cos(qJ(5));
t202 = qJD(5) * t157;
t203 = qJD(5) * t155;
t199 = qJD(1) * qJD(3);
t186 = t158 * t199;
t176 = t152 * t186;
t94 = qJD(3) * t178 - t153 * t176;
t113 = t127 * t152;
t108 = qJD(3) * t113;
t95 = qJD(1) * t108;
t21 = t101 * t202 + t105 * t203 - t155 * t94 + t157 * t95;
t135 = -qJD(5) + t140;
t48 = t101 * t157 + t105 * t155;
t227 = t48 * t135;
t248 = -t21 - t227;
t174 = -t101 * t155 + t105 * t157;
t235 = t174 ^ 2;
t236 = t48 ^ 2;
t247 = t235 - t236;
t234 = t48 * t174;
t246 = -t127 * qJD(3) + t154 * t167;
t173 = t151 * t156 - t153 * t158;
t245 = t140 * t173;
t22 = qJD(5) * t174 - t155 * t95 - t157 * t94;
t228 = t174 * t135;
t244 = -t22 - t228;
t222 = qJ(4) * t152;
t195 = t158 * t222;
t224 = qJ(2) * t156;
t196 = t154 * t224;
t166 = -t195 - t196;
t204 = qJD(4) * t152;
t160 = qJD(3) * t166 - t156 * t204;
t130 = -pkin(2) * t154 - pkin(6) * t152 - pkin(1);
t118 = qJD(1) * t130 + qJD(2);
t200 = qJD(1) * qJD(2);
t187 = t154 * t200;
t205 = qJD(3) * t158;
t214 = t118 * t205 + t158 * t187;
t42 = qJD(1) * t160 + t214;
t223 = qJ(2) * t158;
t139 = t154 * t223;
t209 = qJD(2) * t156;
t188 = t154 * t209;
t165 = -t158 * t204 - t188;
t206 = qJD(3) * t156;
t190 = t118 * t206;
t217 = t152 * t156;
t197 = qJ(4) * t217;
t43 = -t190 + ((-t139 + t197) * qJD(3) + t165) * qJD(1);
t14 = -t151 * t42 + t153 * t43;
t12 = pkin(7) * t95 + t14;
t15 = t151 * t43 + t153 * t42;
t13 = pkin(7) * t94 + t15;
t237 = t105 * pkin(7);
t110 = t158 * t118;
t69 = qJD(1) * t166 + t110;
t58 = -t140 * pkin(3) + t69;
t194 = qJ(2) * t201;
t83 = t118 * t156 + t158 * t194;
t70 = -qJ(4) * t192 + t83;
t62 = t151 * t70;
t30 = t153 * t58 - t62;
t18 = -pkin(4) * t140 - t237 + t30;
t238 = t101 * pkin(7);
t229 = t153 * t70;
t31 = t151 * t58 + t229;
t20 = t31 - t238;
t1 = (qJD(5) * t18 + t13) * t157 + t155 * t12 - t20 * t203;
t119 = pkin(3) * t192 + qJ(2) * t210 + qJD(4);
t68 = pkin(4) * t101 + t119;
t243 = t68 * t48 - t1;
t147 = t152 ^ 2;
t242 = 0.2e1 * t147;
t100 = t130 * t156 + t139;
t6 = t155 * t18 + t157 * t20;
t2 = -qJD(5) * t6 + t12 * t157 - t155 * t13;
t241 = -t174 * t68 + t2;
t240 = t105 ^ 2;
t239 = pkin(3) * t151;
t144 = pkin(3) * t153 + pkin(4);
t115 = t144 * t157 - t155 * t239;
t33 = -t151 * t69 - t229;
t23 = t33 + t238;
t34 = t153 * t69 - t62;
t24 = t34 - t237;
t233 = qJD(5) * t115 - t155 * t23 - t157 * t24;
t116 = t144 * t155 + t157 * t239;
t232 = -qJD(5) * t116 + t155 * t24 - t157 * t23;
t75 = -t127 * t155 - t157 * t173;
t231 = qJD(5) * t75 + t155 * t246 + t157 * t245;
t76 = t127 * t157 - t155 * t173;
t230 = -qJD(5) * t76 - t155 * t245 + t157 * t246;
t208 = qJD(2) * t158;
t213 = t130 * t205 + t154 * t208;
t56 = t160 + t213;
t57 = (-t139 + (-t130 + t222) * t156) * qJD(3) + t165;
t26 = t151 * t57 + t153 * t56;
t125 = t158 * t130;
t74 = -t195 + t125 + (-pkin(3) - t224) * t154;
t84 = -t197 + t100;
t38 = t151 * t74 + t153 * t84;
t177 = t156 * t194;
t82 = t110 - t177;
t226 = t82 * t140;
t225 = t83 * t140;
t221 = t101 * t140;
t220 = t105 * t101;
t219 = t105 * t140;
t159 = qJD(1) ^ 2;
t218 = t147 * t159;
t117 = pkin(3) * t176 + t152 * t200;
t123 = (pkin(3) * t205 + qJD(2)) * t152;
t128 = pkin(3) * t217 + qJ(2) * t152;
t212 = t154 ^ 2 + t147;
t211 = t156 ^ 2 - t158 ^ 2;
t207 = qJD(3) * t152;
t198 = 0.2e1 * t200;
t193 = qJ(2) * t206;
t189 = t119 * t210;
t25 = -t151 * t56 + t153 * t57;
t37 = -t151 * t84 + t153 * t74;
t182 = t212 * t159;
t181 = qJD(1) * t242;
t180 = pkin(3) * t191;
t179 = t158 * t156 * t218;
t61 = -pkin(4) * t94 + t117;
t175 = (-qJD(3) - t140) * t210;
t114 = t173 * t152;
t29 = -pkin(4) * t154 + pkin(7) * t114 + t37;
t32 = -pkin(7) * t113 + t38;
t9 = -t155 * t32 + t157 * t29;
t10 = t155 * t29 + t157 * t32;
t60 = -t113 * t155 - t114 * t157;
t172 = t147 * t156 * t186;
t171 = t212 * t198;
t170 = (t140 + t201) * t207;
t169 = qJ(2) * t205 + t209;
t168 = t154 * t199 + t218;
t164 = -t140 ^ 2 - t218;
t162 = qJD(3) * t101;
t104 = t173 * t207;
t99 = t125 - t196;
t97 = t101 ^ 2;
t81 = -qJD(3) * t100 - t188;
t80 = -t154 * t193 + t213;
t79 = pkin(4) * t105 + t180;
t77 = pkin(4) * t113 + t128;
t71 = -pkin(4) * t104 + t123;
t67 = -t169 * t201 - t190;
t66 = -qJD(3) * t177 + t214;
t59 = t113 * t157 - t114 * t155;
t28 = qJD(5) * t60 - t104 * t157 - t155 * t108;
t27 = -t104 * t155 + t108 * t157 + t113 * t202 - t114 * t203;
t17 = pkin(7) * t104 + t26;
t16 = pkin(7) * t108 + t25;
t5 = -t155 * t20 + t157 * t18;
t4 = -qJD(5) * t10 - t155 * t17 + t157 * t16;
t3 = qJD(5) * t9 + t155 * t16 + t157 * t17;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, qJ(2) * t171, -0.2e1 * t172, t211 * t199 * t242, t156 * t170, 0.2e1 * t172, t158 * t170, 0, -t81 * t140 - t67 * t154 + t169 * t181, t80 * t140 + t66 * t154 + (-t193 + t208) * t181, (-t66 * t156 - t67 * t158 + (t156 * t82 - t158 * t83) * qJD(3) + (-t156 * t80 - t158 * t81 + (-t100 * t158 + t156 * t99) * qJD(3)) * qJD(1)) * t152, qJ(2) * t147 * t198 + t100 * t66 + t67 * t99 + t80 * t83 + t81 * t82, -t105 * t108 + t114 * t95, t101 * t108 + t104 * t105 + t113 * t95 - t114 * t94, t108 * t140 + t154 * t95, -t101 * t104 - t113 * t94, -t104 * t140 - t154 * t94, 0, t101 * t123 - t104 * t119 + t113 * t117 - t128 * t94 - t14 * t154 - t140 * t25, t105 * t123 - t108 * t119 - t114 * t117 - t128 * t95 + t140 * t26 + t15 * t154, -t101 * t26 + t104 * t31 - t105 * t25 + t108 * t30 - t113 * t15 + t114 * t14 + t37 * t95 + t38 * t94, t117 * t128 + t119 * t123 + t14 * t37 + t15 * t38 + t25 * t30 + t26 * t31, -t174 * t27 - t21 * t60, -t174 * t28 + t21 * t59 - t22 * t60 + t27 * t48, t135 * t27 + t154 * t21, t22 * t59 + t28 * t48, t135 * t28 + t154 * t22, 0, -t135 * t4 - t154 * t2 + t22 * t77 + t28 * t68 + t48 * t71 + t59 * t61, t1 * t154 + t135 * t3 + t174 * t71 - t21 * t77 - t27 * t68 + t60 * t61, -t1 * t59 - t10 * t22 - t174 * t4 - t2 * t60 + t21 * t9 + t27 * t5 - t28 * t6 - t3 * t48, t1 * t10 + t2 * t9 + t3 * t6 + t4 * t5 + t61 * t77 + t68 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, -qJ(2) * t182, 0, 0, 0, 0, 0, 0, t164 * t156, t164 * t158, 0, -qJ(2) * t218 + (t67 - t225) * t158 + (t66 + t226) * t156, 0, 0, 0, 0, 0, 0, -t101 * t210 - t140 * t246, -t105 * t210 + t140 * t245, -t101 * t245 - t105 * t246 + t127 * t94 - t173 * t95, t15 * t127 - t14 * t173 + t245 * t31 + t246 * t30 - t189, 0, 0, 0, 0, 0, 0, -t135 * t230 - t210 * t48, t135 * t231 - t174 * t210, -t174 * t230 + t75 * t21 - t76 * t22 - t231 * t48, t1 * t76 + t2 * t75 - t210 * t68 + t230 * t5 + t231 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t211 * t218, t156 * t175, -t179, t158 * t175, 0, -t225 + (-qJD(3) * t118 - t187) * t156 - t168 * t223, t168 * t224 - t214 - t226, 0, 0, t220, -t97 + t240, -t162 - t221, -t220, t94 - t219, 0, -t101 * t180 - t105 * t119 + t140 * t33 + t14, t101 * t119 - t105 * t180 - t140 * t34 - t15, (t31 + t33) * t105 + (-t30 + t34) * t101 + (t151 * t94 + t153 * t95) * pkin(3), -t30 * t33 - t31 * t34 + (t14 * t153 + t15 * t151 - t158 * t189) * pkin(3), t234, t247, t248, -t234, t244, 0, -t135 * t232 - t79 * t48 + t241, t135 * t233 - t174 * t79 + t243, t115 * t21 - t116 * t22 + (-t232 + t6) * t174 + (-t233 - t5) * t48, t1 * t116 + t2 * t115 + t232 * t5 + t233 * t6 - t68 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94 - t219, -t162 + t221, -t97 - t240, t101 * t31 + t105 * t30 + t117, 0, 0, 0, 0, 0, 0, t22 - t228, -t21 + t227, -t235 - t236, t174 * t5 + t48 * t6 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, t247, t248, -t234, t244, 0, -t6 * t135 + t241, -t5 * t135 + t243, 0, 0;];
tauc_reg = t7;
