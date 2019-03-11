% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:07
% EndTime: 2019-03-09 15:30:14
% DurationCPUTime: 2.26s
% Computational Cost: add. (3349->328), mult. (8162->391), div. (0->0), fcn. (5587->6), ass. (0->178)
t144 = sin(qJ(3));
t145 = sin(qJ(2));
t147 = cos(qJ(2));
t234 = cos(qJ(3));
t108 = t144 * t147 + t145 * t234;
t196 = qJD(1) * t108;
t244 = qJD(6) + t196;
t136 = qJD(2) + qJD(3);
t143 = sin(qJ(6));
t146 = cos(qJ(6));
t185 = t234 * t147;
t172 = qJD(1) * t185;
t195 = qJD(1) * t145;
t184 = t144 * t195;
t91 = -t172 + t184;
t68 = t146 * t136 + t143 * t91;
t249 = t244 * t68;
t248 = t146 * t244;
t165 = t143 * t136 - t146 * t91;
t247 = t165 * t244;
t246 = qJD(6) - t244;
t131 = -t147 * pkin(2) - pkin(1);
t237 = -pkin(8) - pkin(7);
t116 = t237 * t147;
t112 = qJD(1) * t116;
t95 = t144 * t112;
t115 = t237 * t145;
t110 = qJD(1) * t115;
t225 = qJD(2) * pkin(2);
t99 = t110 + t225;
t226 = -t234 * t99 - t95;
t245 = qJD(4) + t226;
t191 = qJD(1) * qJD(2);
t243 = -0.2e1 * t191;
t183 = t234 * qJD(3);
t65 = t110 * t234 + t95;
t205 = t183 * pkin(2) + qJD(4) - t65;
t194 = qJD(3) * t144;
t98 = t234 * t112;
t64 = t144 * t110 - t98;
t170 = pkin(2) * t194 - t64;
t90 = t196 ^ 2;
t242 = -t136 ^ 2 - t90;
t204 = t144 * t145;
t168 = t136 * t204;
t198 = t136 * t172;
t55 = qJD(1) * t168 - t198;
t218 = t146 * t55;
t219 = t143 * t244;
t157 = t219 * t244 + t218;
t67 = t136 * t108;
t56 = t67 * qJD(1);
t27 = -qJD(6) * t165 + t143 * t56;
t186 = qJD(2) * t237;
t173 = qJD(1) * t186;
t100 = t145 * t173;
t101 = t147 * t173;
t24 = t144 * t100 - t234 * t101 - t112 * t183 + t99 * t194;
t162 = t55 * qJ(5) + t24;
t13 = -qJD(5) * t196 + t162;
t111 = t145 * t186;
t113 = t147 * t186;
t72 = t144 * t115 - t116 * t234;
t32 = qJD(3) * t72 + t144 * t111 - t113 * t234;
t66 = -qJD(2) * t185 - t147 * t183 + t168;
t16 = t66 * qJ(5) - t108 * qJD(5) + t32;
t107 = -t185 + t204;
t133 = t136 * qJD(4);
t160 = -t234 * t100 - t144 * t101 - t112 * t194 - t99 * t183;
t22 = t133 - t160;
t12 = -t56 * qJ(5) - t91 * qJD(5) - t22;
t71 = -t115 * t234 - t144 * t116;
t48 = -t108 * qJ(5) + t71;
t169 = -t12 * t107 + t48 * t55;
t114 = t131 * qJD(1);
t45 = t91 * pkin(3) - qJ(4) * t196 + t114;
t167 = qJD(5) - t45;
t238 = -pkin(4) - pkin(9);
t18 = pkin(5) * t196 + t238 * t91 + t167;
t174 = -t107 * pkin(3) + t108 * qJ(4) - t131;
t23 = t108 * pkin(5) + t107 * t238 + t174;
t134 = t136 * qJ(4);
t210 = t91 * qJ(5);
t60 = t144 * t99 - t98;
t39 = t210 + t60;
t34 = -t134 - t39;
t30 = t136 * pkin(5) - t34;
t240 = -(qJD(6) * t18 + t13) * t108 - (qJD(6) * t23 + t16) * t244 + t30 * t67 + t169;
t148 = -pkin(3) - pkin(4);
t138 = -pkin(3) + t238;
t208 = t196 * qJ(5);
t164 = -t208 + t245;
t25 = t136 * t138 + t164;
t6 = t143 * t18 + t146 * t25;
t236 = t6 * t91;
t235 = t196 * pkin(4);
t233 = t23 * t55;
t232 = t45 * t91;
t231 = t45 * t196;
t230 = t244 * t91;
t229 = t91 * t68;
t228 = t91 * t165;
t227 = t196 * t91;
t57 = pkin(3) * t196 + t91 * qJ(4);
t224 = t114 * t91;
t223 = t114 * t196;
t128 = t144 * pkin(2) + qJ(4);
t222 = t128 * t56;
t221 = t143 * t55;
t217 = t146 * t165;
t192 = qJD(6) * t146;
t193 = qJD(6) * t143;
t26 = -t136 * t192 + t146 * t56 - t193 * t91;
t215 = t26 * t143;
t214 = t30 * t107;
t31 = t234 * t111 + t144 * t113 + t115 * t183 + t116 * t194;
t213 = t31 * t136;
t212 = t32 * t136;
t211 = t56 * qJ(4);
t209 = t91 * t136;
t206 = -t208 + t205;
t150 = qJD(1) ^ 2;
t203 = t147 * t150;
t149 = qJD(2) ^ 2;
t202 = t149 * t145;
t201 = t149 * t147;
t197 = t145 ^ 2 - t147 ^ 2;
t189 = t145 * t225;
t188 = pkin(2) * t195;
t166 = t143 * t25 - t146 * t18;
t187 = -t12 * t146 - t166 * t91;
t182 = t145 * t191;
t181 = t244 * t30;
t177 = pkin(1) * t243;
t130 = -pkin(2) * t234 - pkin(3);
t127 = -pkin(4) + t130;
t121 = -pkin(9) + t127;
t20 = -t91 * pkin(5) + t196 * t238 - t57;
t176 = qJD(6) * t121 - t188 + t20;
t175 = qJD(6) * t138 + t20;
t171 = -t210 + t170;
t47 = t188 + t57;
t163 = -t24 - t231;
t161 = t60 * t136 - t24;
t21 = t67 * pkin(3) + t66 * qJ(4) - t108 * qJD(4) + t189;
t17 = pkin(2) * t182 + t56 * pkin(3) + t55 * qJ(4) - qJD(4) * t196;
t158 = -t244 * t248 + t221;
t156 = -t136 * t184 + t198;
t155 = t138 * t55 + t244 * t39 - t181;
t154 = -t136 * t226 + t160;
t10 = -t56 * pkin(4) - t17;
t29 = -t91 * pkin(4) + t167;
t153 = (-qJD(5) - t29) * t196 + t162;
t152 = t29 * t91 - t12;
t151 = t121 * t55 - t171 * t244 - t181;
t142 = qJ(4) + pkin(5);
t125 = pkin(5) + t128;
t89 = t91 ^ 2;
t50 = t134 + t60;
t49 = t107 * qJ(5) + t72;
t46 = -t136 * pkin(3) + t245;
t44 = -t89 + t90;
t41 = -t107 * pkin(4) + t174;
t40 = t55 * t108;
t37 = -t57 - t235;
t35 = t156 + t209;
t33 = -t47 - t235;
t28 = t136 * t148 + t164;
t15 = -t67 * qJ(5) - t107 * qJD(5) - t31;
t14 = -t67 * pkin(4) - t21;
t9 = t158 - t228;
t8 = t157 - t229;
t7 = t165 * t248 - t215;
t4 = -t66 * pkin(5) + t238 * t67 - t21;
t3 = -t55 * pkin(5) + t238 * t56 - t17;
t2 = t146 * t3;
t1 = (-t26 + t249) * t146 + (t27 - t247) * t143;
t5 = [0, 0, 0, 0.2e1 * t147 * t182, t197 * t243, t201, -t202, 0, -pkin(7) * t201 + t145 * t177, pkin(7) * t202 + t147 * t177, -t196 * t66 - t40, t55 * t107 - t108 * t56 - t196 * t67 + t66 * t91, -t66 * t136, -t67 * t136, 0, t114 * t67 + t131 * t56 - t212 + (qJD(1) * t107 + t91) * t189, -t114 * t66 - t131 * t55 + 0.2e1 * t196 * t189 - t213, t17 * t107 - t174 * t56 + t21 * t91 + t45 * t67 - t212, -t22 * t107 + t24 * t108 + t196 * t32 - t31 * t91 - t46 * t66 - t50 * t67 - t71 * t55 - t72 * t56, -t17 * t108 - t174 * t55 - t196 * t21 + t45 * t66 + t213, -t17 * t174 + t45 * t21 + t22 * t72 + t24 * t71 + t50 * t31 + t46 * t32, t10 * t108 - t15 * t136 + t14 * t196 - t29 * t66 - t41 * t55, t10 * t107 + t16 * t136 + t14 * t91 + t29 * t67 + t41 * t56, -t13 * t108 - t15 * t91 - t16 * t196 + t28 * t66 - t34 * t67 + t49 * t56 + t169, t10 * t41 - t12 * t49 + t13 * t48 + t29 * t14 + t34 * t15 + t28 * t16, -t67 * t217 + (t26 * t146 + t165 * t193) * t107 (t143 * t165 - t146 * t68) * t67 + (-t215 - t146 * t27 + (t143 * t68 + t217) * qJD(6)) * t107, t67 * t248 + t26 * t108 + t165 * t66 + (-t193 * t244 - t218) * t107, -t67 * t219 - t27 * t108 + t68 * t66 + (-t192 * t244 + t221) * t107, -t244 * t66 - t40, t2 * t108 - t15 * t68 + t49 * t27 + t166 * t66 + (-t233 + t4 * t244 + (-t25 * t108 - t244 * t48 + t214) * qJD(6)) * t146 + t240 * t143, t15 * t165 + t49 * t26 + t6 * t66 + (-(-qJD(6) * t48 + t4) * t244 + t233 - (-qJD(6) * t25 + t3) * t108 - qJD(6) * t214) * t143 + t240 * t146; 0, 0, 0, -t145 * t203, t197 * t150, 0, 0, 0, t150 * pkin(1) * t145, pkin(1) * t203, t227, t44, t35, 0, 0, -t223 + t64 * t136 + (-t136 * t194 - t195 * t91) * pkin(2) - t24, t224 + t65 * t136 + (-t136 * t183 - t195 * t196) * pkin(2) + t160, -t136 * t170 - t47 * t91 + t163, -t222 - t130 * t55 + (t170 + t50) * t196 + (t46 - t205) * t91, t136 * t205 + t196 * t47 + t22 - t232, t22 * t128 + t24 * t130 + t170 * t46 + t205 * t50 - t45 * t47, t136 * t206 - t196 * t33 + t152, t136 * t171 - t33 * t91 + t153, t127 * t55 + t222 + (-t171 + t34) * t196 + (-t28 + t206) * t91, -t12 * t128 + t13 * t127 + t171 * t28 - t206 * t34 - t29 * t33, t7, t1, t9, t8, t230, t125 * t27 + t143 * t151 - t176 * t248 + t206 * t68 + t187, t125 * t26 - t236 - t206 * t165 + (t176 * t244 + t12) * t143 + t151 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t44, t35, 0, 0, t161 - t223, t154 + t224, -t57 * t91 + t161 - t231, pkin(3) * t55 - t211 + (t50 - t60) * t196 + (t46 - t245) * t91, t196 * t57 + 0.2e1 * t133 - t154 - t232, -t24 * pkin(3) + t22 * qJ(4) + t245 * t50 - t45 * t57 - t46 * t60, t136 * t164 - t196 * t37 + t152, -t39 * t136 - t37 * t91 + t153, t211 + t148 * t55 + (t34 + t39) * t196 + (-t28 + t164) * t91, -t12 * qJ(4) + t13 * t148 - t164 * t34 - t28 * t39 - t29 * t37, t7, t1, t9, t8, t230, t142 * t27 + t143 * t155 + t164 * t68 - t175 * t248 + t187, t142 * t26 - t236 - t164 * t165 + (t175 * t244 + t12) * t143 + t155 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t35, t242, -t50 * t136 - t163, t242, -t227, t55 - t209, t34 * t136 + t153, 0, 0, 0, 0, 0, -t136 * t68 + t158, t136 * t165 + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156 - t209, t136 * t196 + t56, -t90 - t89, t196 * t28 + t34 * t91 + t10, 0, 0, 0, 0, 0, -t157 - t229, t158 + t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165 * t68, t165 ^ 2 - t68 ^ 2, t26 + t249, -t27 - t247, -t55, -t143 * t13 + t30 * t165 - t246 * t6 + t2, -t146 * t13 - t143 * t3 + t246 * t166 + t30 * t68;];
tauc_reg  = t5;
