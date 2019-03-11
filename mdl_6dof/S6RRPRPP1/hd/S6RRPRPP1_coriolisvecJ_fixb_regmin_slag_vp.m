% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:47:52
% EndTime: 2019-03-09 09:48:00
% DurationCPUTime: 2.98s
% Computational Cost: add. (7802->380), mult. (19785->508), div. (0->0), fcn. (14409->8), ass. (0->192)
t182 = sin(pkin(9));
t184 = cos(pkin(9));
t186 = sin(qJ(2));
t188 = cos(qJ(2));
t160 = t182 * t188 + t184 * t186;
t181 = sin(pkin(10));
t183 = cos(pkin(10));
t185 = sin(qJ(4));
t187 = cos(qJ(4));
t261 = -t181 * t185 + t183 * t187;
t99 = t261 * t160;
t224 = qJD(4) * t185;
t225 = qJD(1) * t186;
t237 = t184 * t188;
t146 = qJD(1) * t237 - t182 * t225;
t244 = t146 * t185;
t268 = t224 - t244;
t141 = qJD(4) - t146;
t149 = t160 * qJD(1);
t222 = t187 * qJD(2);
t119 = t185 * t149 - t222;
t121 = t185 * qJD(2) + t187 * t149;
t66 = t183 * t119 + t181 * t121;
t267 = t141 * t66;
t159 = t181 * t187 + t183 * t185;
t145 = t159 * qJD(4);
t249 = -t159 * t146 + t145;
t223 = qJD(4) * t187;
t248 = -t261 * t146 - t181 * t224 + t183 * t223;
t215 = t160 * t223;
t158 = t182 * t186 - t237;
t151 = t158 * qJD(2);
t234 = t185 * t151;
t266 = t215 - t234;
t201 = -t181 * t119 + t183 * t121;
t265 = t201 ^ 2;
t221 = qJD(1) * qJD(2);
t264 = -0.2e1 * t221;
t174 = t182 * pkin(2) + pkin(8);
t229 = qJ(5) + t174;
t209 = qJD(4) * t229;
t126 = t187 * qJD(5) - t185 * t209;
t193 = -t185 * qJD(5) - t187 * t209;
t258 = -qJ(3) - pkin(7);
t166 = t258 * t188;
t163 = qJD(1) * t166;
t152 = t182 * t163;
t165 = t258 * t186;
t162 = qJD(1) * t165;
t112 = t184 * t162 + t152;
t94 = pkin(2) * t225 + t149 * pkin(3) - t146 * pkin(8);
t90 = t187 * t94;
t46 = -t187 * t146 * qJ(5) + t149 * pkin(4) - t185 * t112 + t90;
t250 = t187 * t112 + t185 * t94;
t55 = -qJ(5) * t244 + t250;
t254 = (t193 - t46) * t183 + (-t126 + t55) * t181;
t253 = qJD(2) * pkin(2);
t155 = t162 + t253;
t109 = t184 * t155 + t152;
t102 = -qJD(2) * pkin(3) - t109;
t64 = t119 * pkin(4) + qJD(5) + t102;
t35 = t66 * pkin(5) - qJ(6) * t201 + t64;
t263 = t35 * t201;
t238 = t184 * t163;
t111 = t182 * t162 - t238;
t262 = pkin(4) * t268 - t111;
t213 = t188 * t221;
t214 = t186 * t221;
t137 = -t182 * t214 + t184 * t213;
t81 = qJD(4) * t121 + t185 * t137;
t260 = t141 ^ 2;
t110 = t182 * t155 - t238;
t103 = qJD(2) * pkin(8) + t110;
t216 = -t188 * pkin(2) - pkin(1);
t203 = t216 * qJD(1);
t164 = qJD(3) + t203;
t84 = -t146 * pkin(3) - t149 * pkin(8) + t164;
t58 = t187 * t103 + t185 * t84;
t48 = -t119 * qJ(5) + t58;
t43 = t183 * t48;
t57 = -t185 * t103 + t187 * t84;
t47 = -t121 * qJ(5) + t57;
t23 = t181 * t47 + t43;
t259 = t23 * t201;
t147 = t160 * qJD(2);
t136 = qJD(1) * t147;
t170 = pkin(2) * t214;
t83 = t136 * pkin(3) - t137 * pkin(8) + t170;
t76 = t187 * t83;
t212 = qJD(2) * t258;
t143 = t188 * qJD(3) + t186 * t212;
t129 = t143 * qJD(1);
t144 = -t186 * qJD(3) + t188 * t212;
t130 = t144 * qJD(1);
t79 = t184 * t129 + t182 * t130;
t191 = -qJD(4) * t58 - t185 * t79 + t76;
t80 = qJD(4) * t222 + t187 * t137 - t149 * t224;
t13 = t136 * pkin(4) - t80 * qJ(5) - t121 * qJD(5) + t191;
t195 = -t103 * t224 + t185 * t83 + t187 * t79 + t84 * t223;
t18 = -t81 * qJ(5) - t119 * qJD(5) + t195;
t3 = t183 * t13 - t181 * t18;
t4 = t181 * t13 + t183 * t18;
t28 = t181 * t46 + t183 * t55;
t21 = t149 * qJ(6) + t28;
t74 = t183 * t126 + t181 * t193;
t257 = -t21 + t74;
t256 = -t149 * pkin(5) + t254;
t108 = t158 * pkin(3) - t160 * pkin(8) + t216;
t117 = t182 * t165 - t184 * t166;
t113 = t187 * t117;
t200 = qJ(5) * t151 - qJD(5) * t160;
t218 = t186 * t253;
t95 = t147 * pkin(3) + t151 * pkin(8) + t218;
t91 = t187 * t95;
t93 = t184 * t143 + t182 * t144;
t26 = t147 * pkin(4) - t185 * t93 + t91 + t200 * t187 + (-t113 + (qJ(5) * t160 - t108) * t185) * qJD(4);
t220 = t108 * t223 + t185 * t95 + t187 * t93;
t30 = -qJ(5) * t215 + (-qJD(4) * t117 + t200) * t185 + t220;
t8 = t181 * t26 + t183 * t30;
t255 = -t249 * pkin(5) + t248 * qJ(6) + t159 * qJD(6) - t262;
t40 = t141 * pkin(4) + t47;
t20 = t181 * t40 + t43;
t101 = t187 * t108;
t241 = t160 * t187;
t54 = t158 * pkin(4) - qJ(5) * t241 - t185 * t117 + t101;
t227 = t185 * t108 + t113;
t242 = t160 * t185;
t60 = -qJ(5) * t242 + t227;
t34 = t181 * t54 + t183 * t60;
t252 = t181 * t48;
t251 = t80 * t185;
t247 = t119 * t141;
t246 = t121 * t141;
t245 = t121 * t149;
t243 = t149 * t119;
t236 = t185 * t136;
t127 = t187 * t136;
t233 = t187 * t151;
t190 = qJD(1) ^ 2;
t232 = t188 * t190;
t189 = qJD(2) ^ 2;
t231 = t189 * t186;
t230 = t189 * t188;
t24 = t183 * t47 - t252;
t228 = qJD(6) - t24;
t78 = t182 * t129 - t184 * t130;
t226 = t186 ^ 2 - t188 ^ 2;
t219 = t136 * qJ(6) + t4;
t176 = -t184 * pkin(2) - pkin(3);
t50 = t181 * t80 + t183 * t81;
t211 = t229 * t185;
t210 = pkin(1) * t264;
t92 = t182 * t143 - t184 * t144;
t116 = -t184 * t165 - t182 * t166;
t208 = t141 * t187;
t2 = -t136 * pkin(5) - t3;
t56 = t81 * pkin(4) + t78;
t156 = t229 * t187;
t104 = t181 * t156 + t183 * t211;
t105 = t183 * t156 - t181 * t211;
t51 = -t181 * t81 + t183 * t80;
t207 = t104 * t51 - t105 * t50 - t74 * t66;
t206 = -t66 ^ 2 - t265;
t204 = pkin(4) * t242 + t116;
t7 = -t181 * t30 + t183 * t26;
t19 = t183 * t40 - t252;
t33 = -t181 * t60 + t183 * t54;
t202 = -t117 * t136 + t78 * t160;
t199 = -t187 * pkin(4) + t176;
t198 = t266 * pkin(4) + t92;
t197 = -t141 * t268 + t127;
t196 = -t160 * t224 - t233;
t194 = t102 * t141 - t174 * t136;
t192 = -t159 * t50 + t201 * t249 - t248 * t66 - t261 * t51;
t9 = t50 * pkin(5) - t51 * qJ(6) - qJD(6) * t201 + t56;
t175 = -t183 * pkin(4) - pkin(5);
t171 = t181 * pkin(4) + qJ(6);
t98 = t159 * t160;
t96 = -pkin(5) * t261 - t159 * qJ(6) + t199;
t62 = t145 * t160 - t181 * t234 + t183 * t233;
t61 = -qJD(4) * t99 + t151 * t159;
t41 = t98 * pkin(5) - t99 * qJ(6) + t204;
t37 = t121 * pkin(4) + pkin(5) * t201 + qJ(6) * t66;
t32 = -t158 * pkin(5) - t33;
t31 = t158 * qJ(6) + t34;
t15 = t141 * qJ(6) + t20;
t14 = -t141 * pkin(5) + qJD(6) - t19;
t12 = -t61 * pkin(5) + t62 * qJ(6) - t99 * qJD(6) + t198;
t6 = -t147 * pkin(5) - t7;
t5 = t147 * qJ(6) + t158 * qJD(6) + t8;
t1 = t141 * qJD(6) + t219;
t10 = [0, 0, 0, 0.2e1 * t186 * t213, t226 * t264, t230, -t231, 0, -pkin(7) * t230 + t186 * t210, pkin(7) * t231 + t188 * t210, t109 * t151 - t110 * t147 + t116 * t137 + t93 * t146 + t92 * t149 - t79 * t158 + t202, -t109 * t92 + t110 * t93 + t78 * t116 + t79 * t117 + (t164 + t203) * t218, t121 * t196 + t80 * t241 -(-t119 * t187 - t121 * t185) * t151 + (-t251 - t187 * t81 + (t119 * t185 - t121 * t187) * qJD(4)) * t160, t121 * t147 + t160 * t127 + t141 * t196 + t80 * t158, -t119 * t147 - t266 * t141 - t81 * t158 - t160 * t236, t136 * t158 + t141 * t147 (-t117 * t223 + t91) * t141 + t101 * t136 + (-t103 * t223 + t76) * t158 + t57 * t147 + t92 * t119 + t116 * t81 + t102 * t215 + ((-qJD(4) * t108 - t93) * t141 + (-qJD(4) * t84 - t79) * t158 - t102 * t151 + t202) * t185 -(-t117 * t224 + t220) * t141 - t227 * t136 - t195 * t158 - t58 * t147 + t92 * t121 + t116 * t80 + t78 * t241 + t196 * t102, t19 * t62 + t20 * t61 - t201 * t7 - t3 * t99 - t33 * t51 - t34 * t50 - t4 * t98 - t8 * t66, t19 * t7 + t198 * t64 + t20 * t8 + t204 * t56 + t3 * t33 + t4 * t34, t12 * t66 - t32 * t136 - t14 * t147 - t6 * t141 - t2 * t158 - t35 * t61 + t41 * t50 + t9 * t98, -t1 * t98 - t14 * t62 + t15 * t61 + t2 * t99 + t201 * t6 - t31 * t50 + t32 * t51 - t5 * t66, t1 * t158 - t12 * t201 + t31 * t136 + t5 * t141 + t15 * t147 + t35 * t62 - t41 * t51 - t9 * t99, t1 * t31 + t35 * t12 + t14 * t6 + t15 * t5 + t2 * t32 + t9 * t41; 0, 0, 0, -t186 * t232, t226 * t190, 0, 0, 0, t190 * pkin(1) * t186, pkin(1) * t232 (t110 - t111) * t149 + (t109 - t112) * t146 + (-t136 * t182 - t137 * t184) * pkin(2), t109 * t111 - t110 * t112 + (-t164 * t225 + t182 * t79 - t184 * t78) * pkin(2), t121 * t208 + t251 (t80 - t247) * t187 + (-t81 - t246) * t185, t141 * t208 + t236 - t245, t197 + t243, -t141 * t149, -t111 * t119 - t57 * t149 + t176 * t81 - t78 * t187 + (-t174 * t223 - t90) * t141 + (t112 * t141 + t194) * t185, -t111 * t121 + t58 * t149 + t176 * t80 + t78 * t185 + (t174 * t224 + t250) * t141 + t194 * t187, -t3 * t159 - t248 * t19 - t249 * t20 - t201 * t254 + t261 * t4 + t28 * t66 + t207, t4 * t105 - t3 * t104 + t56 * t199 + t262 * t64 + (t74 - t28) * t20 + t254 * t19, -t104 * t136 + t14 * t149 + t256 * t141 + t249 * t35 - t255 * t66 - t261 * t9 + t96 * t50, t1 * t261 + t248 * t14 - t249 * t15 + t2 * t159 - t201 * t256 + t21 * t66 + t207, t105 * t136 + t257 * t141 - t15 * t149 - t9 * t159 + t201 * t255 - t248 * t35 - t96 * t51, t1 * t105 + t2 * t104 - t256 * t14 + t257 * t15 - t255 * t35 + t9 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146 ^ 2 - t149 ^ 2, t109 * t149 - t110 * t146 + t170, 0, 0, 0, 0, 0, t197 - t243, -t260 * t187 - t236 - t245, t192, -t64 * t149 + t4 * t159 - t249 * t19 + t248 * t20 + t261 * t3, t136 * t261 - t249 * t141 - t149 * t66, t192, t159 * t136 + t248 * t141 + t149 * t201, t1 * t159 + t14 * t249 - t35 * t149 + t15 * t248 - t2 * t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t119, -t119 ^ 2 + t121 ^ 2, t80 + t247, t246 - t81, t136, -t102 * t121 + t58 * t141 + t191, t102 * t119 + t57 * t141 - t195, t20 * t201 - t259 + (-t181 * t50 - t183 * t51) * pkin(4) + (-t19 + t24) * t66, t19 * t23 - t20 * t24 + (-t121 * t64 + t181 * t4 + t183 * t3) * pkin(4), t23 * t141 - t263 - t37 * t66 + (pkin(5) - t175) * t136 + t3, t15 * t201 - t171 * t50 + t175 * t51 - t259 + (t14 - t228) * t66, t171 * t136 - t35 * t66 + t37 * t201 + (0.2e1 * qJD(6) - t24) * t141 + t219, t1 * t171 - t14 * t23 + t15 * t228 + t2 * t175 - t35 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, t19 * t201 + t20 * t66 + t56, t141 * t201 + t50, t206, -t51 + t267, -t14 * t201 + t15 * t66 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t149 + t201 * t66, t51 + t267, -t260 - t265, -t15 * t141 + t2 + t263;];
tauc_reg  = t10;
