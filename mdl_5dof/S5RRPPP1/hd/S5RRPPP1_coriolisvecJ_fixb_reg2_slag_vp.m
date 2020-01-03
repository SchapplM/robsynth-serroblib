% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:31
% EndTime: 2019-12-31 19:24:40
% DurationCPUTime: 2.75s
% Computational Cost: add. (2982->363), mult. (9489->505), div. (0->0), fcn. (6891->6), ass. (0->180)
t151 = cos(qJ(2));
t149 = cos(pkin(5));
t201 = t149 * qJD(2);
t148 = cos(pkin(8));
t207 = qJD(2) * t148;
t146 = sin(pkin(8));
t150 = sin(qJ(2));
t219 = t150 * t146;
t100 = t151 * t207 - t201 * t219;
t147 = sin(pkin(5));
t209 = qJD(1) * t151;
t122 = t147 * t209 - t201;
t206 = qJD(2) * t150;
t218 = t150 * t148;
t220 = t149 * t151;
t110 = t146 * t220 + t218;
t211 = qJD(1) * t110;
t208 = qJD(2) * t147;
t74 = t146 * t208 + t211;
t210 = qJD(1) * t150;
t196 = t146 * t210;
t185 = t149 * t196;
t200 = qJD(1) * qJD(2);
t192 = t151 * t200;
t88 = -qJD(2) * t185 + t148 * t192;
t12 = -t100 * t122 + ((t74 + t211) * t206 - t151 * t88) * t147;
t109 = -t148 * t220 + t219;
t195 = t149 * t209;
t73 = -t147 * t207 - t148 * t195 + t196;
t166 = t146 * t151 + t149 * t218;
t99 = t166 * qJD(2);
t87 = qJD(1) * t99;
t11 = -t122 * t99 + ((qJD(1) * t109 + t73) * t206 - t151 * t87) * t147;
t248 = -0.2e1 * t200;
t197 = t147 * t210;
t203 = qJD(3) * t148;
t222 = t147 * t151;
t170 = pkin(2) * t150 - qJ(3) * t222;
t116 = t170 * qJD(1);
t190 = qJ(3) * t149 + pkin(7);
t174 = qJD(1) * t190;
t117 = t150 * t174;
t118 = t151 * t174;
t224 = t146 * t149;
t225 = t146 * t147;
t45 = t116 * t225 - t148 * t117 - t118 * t224;
t229 = qJ(4) * t197 - t149 * qJD(4) - t147 * t203 + t45;
t119 = t122 ^ 2;
t247 = -t74 ^ 2 - t119;
t125 = -t150 * t147 * qJ(3) - t151 * pkin(2) - pkin(1);
t126 = t190 * t150;
t246 = t148 * (t125 * t147 - t126 * t149);
t143 = t147 ^ 2;
t245 = t149 ^ 2 + t143;
t202 = qJD(3) * t150;
t84 = qJD(2) * t170 - t147 * t202;
t173 = qJD(2) * t190;
t85 = qJD(3) * t220 - t150 * t173;
t86 = -t149 * t202 - t151 * t173;
t31 = t148 * t85 + t86 * t224 + t84 * t225;
t21 = -t147 * (qJ(4) * t206 - qJD(4) * t151) - t31;
t105 = qJD(2) * pkin(2) - t117;
t106 = t125 * qJD(1);
t124 = t195 + t208;
t92 = pkin(7) * t209 + qJ(3) * t124;
t38 = -t146 * t92 + (t105 * t149 + t106 * t147) * t148;
t97 = t166 * qJD(1);
t98 = t148 * t209 - t185;
t244 = -t73 * t98 - t97 * t74 + (t146 * t87 - t148 * t88) * t147;
t243 = (t143 * t207 + t147 * t73) * t210 - t122 * t97 - t149 * t87;
t16 = (qJD(2) * t143 * t146 - t147 * t74) * t210 + t98 * t122 + t88 * t149;
t239 = t74 * t73;
t238 = pkin(3) + qJ(5);
t71 = t84 * qJD(1);
t236 = t148 * t71;
t235 = t148 * t84;
t52 = t245 * t148 * t210 + t146 * t124;
t232 = t52 * t122;
t54 = t149 * t116 + t147 * t118;
t169 = -t98 * qJ(4) + t54;
t231 = -t238 * t97 + (-qJD(4) * t146 - qJD(5) * t148) * t147 - t169;
t102 = t146 * t117;
t221 = t148 * t149;
t180 = t118 * t221 - t102;
t194 = t238 * t150;
t204 = qJD(3) * t147;
t226 = t116 * t148;
t230 = -t149 * qJD(5) + t146 * t204 - t98 * pkin(4) - (-qJD(1) * t194 - t226) * t147 - t180;
t228 = t97 * pkin(4) - t229;
t223 = t147 * t148;
t153 = qJD(1) ^ 2;
t217 = t151 * t153;
t152 = qJD(2) ^ 2;
t216 = t152 * t150;
t215 = t152 * t151;
t214 = t74 * qJD(4);
t127 = t190 * t151;
t113 = t146 * t127;
t213 = pkin(3) * t222 + t113;
t112 = pkin(2) * t224 + qJ(3) * t223;
t212 = t150 ^ 2 - t151 ^ 2;
t205 = qJD(3) * t122;
t61 = qJD(1) * t85 + qJD(2) * t204;
t72 = t86 * qJD(1);
t19 = t148 * t61 + t72 * t224 + t71 * t225;
t39 = t105 * t224 + t106 * t225 + t148 * t92;
t49 = t125 * t225 - t126 * t224 + t148 * t127;
t199 = t150 * t217;
t198 = -pkin(2) * t148 - pkin(3);
t193 = t150 * t200;
t191 = -qJ(4) * t146 - pkin(2);
t136 = t147 * t193;
t189 = -t136 + t239;
t47 = -t147 * t72 + t149 * t71;
t50 = -t147 * t86 + t149 * t84;
t188 = pkin(1) * t248;
t58 = t149 * t125 + t147 * t126;
t184 = t150 * t192;
t55 = t146 * t61;
t183 = -t72 * t221 + t55;
t76 = t146 * t85;
t182 = -t86 * t221 + t76;
t181 = qJD(2) * t194;
t9 = -qJ(4) * t136 + t122 * qJD(4) - t19;
t53 = t148 * t124 - t245 * t196;
t179 = -t52 * t74 + t53 * t73;
t93 = -t149 * qJ(4) - t112;
t51 = -t147 * t105 + t149 * t106 + qJD(3);
t32 = t87 * t109 + t73 * t99;
t33 = t74 * t100 + t88 * t110;
t27 = t122 * qJ(4) - t39;
t37 = t88 * t225 - t74 * t98;
t36 = -t87 * t223 - t73 * t97;
t165 = -t110 * qJ(4) + t58;
t43 = qJ(4) * t222 - t49;
t3 = -t87 * pkin(4) - t9;
t163 = -t74 * qJ(4) + t51;
t162 = -t53 * t122 + t88;
t161 = -t122 * t73 + t88;
t160 = -t100 * qJ(4) - t110 * qJD(4) + t50;
t159 = t166 * t200;
t158 = qJD(4) - t38;
t157 = t73 * t100 + t109 * t88 + t87 * t110 + t99 * t74;
t5 = t87 * pkin(3) - t88 * qJ(4) - t214 + t47;
t155 = t87 - t232;
t1 = t87 * qJ(5) + t73 * qJD(5) + t5;
t14 = (-pkin(3) * t193 - t236) * t147 + t183;
t154 = t88 * pkin(4) + (-qJD(1) * t181 - t236) * t147 + t183;
t138 = qJ(3) * t225;
t111 = pkin(2) * t221 - t138;
t95 = (-pkin(3) * t148 + t191) * t147;
t94 = t149 * t198 + t138;
t69 = (-t238 * t148 + t191) * t147;
t68 = pkin(4) * t223 - t93;
t63 = pkin(4) * t225 + t138 + (-qJ(5) + t198) * t149;
t62 = (-t122 * t147 - t143 * t209) * t206;
t57 = (t122 + t201) * t197;
t48 = -t113 + t246;
t46 = t213 - t246;
t44 = t102 + (t116 * t147 - t118 * t149) * t148;
t42 = t109 * pkin(3) + t165;
t41 = (-pkin(3) * t210 - t226) * t147 + t180;
t35 = t97 * pkin(3) + t169;
t34 = -t109 * pkin(4) - t43;
t30 = -t76 + (t147 * t84 + t149 * t86) * t148;
t29 = t238 * t109 + t165;
t28 = t126 * t221 + t110 * pkin(4) + (qJ(5) * t151 - t125 * t148) * t147 + t213;
t26 = t122 * pkin(3) + t158;
t24 = (-pkin(3) * t206 - t235) * t147 + t182;
t20 = t73 * pkin(3) + t163;
t18 = -t55 + (t147 * t71 + t149 * t72) * t148;
t15 = t99 * pkin(3) + t160;
t13 = -t73 * pkin(4) + qJD(5) - t27;
t10 = -t99 * pkin(4) - t21;
t8 = t100 * pkin(4) + (qJD(5) * t151 - t181 - t235) * t147 + t182;
t7 = t238 * t73 + t163;
t6 = t74 * pkin(4) + t238 * t122 + t158;
t4 = t109 * qJD(5) + t238 * t99 + t160;
t2 = t122 * qJD(5) + t154;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t184, t212 * t248, t215, -0.2e1 * t184, -t216, 0, -pkin(7) * t215 + t150 * t188, pkin(7) * t216 + t151 * t188, 0, 0, t33, -t157, t12, t32, -t11, t62, t47 * t109 - t30 * t122 + t50 * t73 + t51 * t99 + t58 * t87 + (-t151 * t18 + (qJD(1) * t48 + t38) * t206) * t147, t51 * t100 + t47 * t110 + t31 * t122 + t50 * t74 + t58 * t88 + (t151 * t19 + (-qJD(1) * t49 - t39) * t206) * t147, -t38 * t100 - t19 * t109 - t18 * t110 - t30 * t74 - t31 * t73 - t39 * t99 - t48 * t88 - t49 * t87, t18 * t48 + t19 * t49 + t38 * t30 + t39 * t31 + t47 * t58 + t51 * t50, t62, -t12, t11, t33, -t157, t32, t26 * t100 + t9 * t109 + t14 * t110 + t21 * t73 + t24 * t74 + t27 * t99 + t43 * t87 + t46 * t88, -t5 * t109 - t24 * t122 - t15 * t73 - t20 * t99 - t42 * t87 + (-t14 * t151 + (qJD(1) * t46 + t26) * t206) * t147, -t20 * t100 - t5 * t110 + t21 * t122 - t15 * t74 - t42 * t88 + (t151 * t9 + (-qJD(1) * t43 - t27) * t206) * t147, t14 * t46 + t20 * t15 + t27 * t21 + t26 * t24 + t5 * t42 + t9 * t43, t62, t11, t12, t32, t157, t33, -t10 * t73 + t6 * t100 - t3 * t109 + t2 * t110 - t13 * t99 + t28 * t88 - t34 * t87 + t8 * t74, -t1 * t110 - t10 * t122 - t7 * t100 - t29 * t88 - t4 * t74 + (-t151 * t3 + (qJD(1) * t34 + t13) * t206) * t147, t1 * t109 + t8 * t122 + t29 * t87 + t4 * t73 + t7 * t99 + (t151 * t2 + (-qJD(1) * t28 - t6) * t206) * t147, t1 * t29 + t13 * t10 + t2 * t28 + t3 * t34 + t7 * t4 + t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, t212 * t153, 0, t199, 0, 0, t153 * pkin(1) * t150, pkin(1) * t217, 0, 0, t37, -t244, t16, t36, t243, t57, t44 * t122 + t18 * t149 - t51 * t97 - t54 * t73 + (t146 * t205 - pkin(2) * t87 - t148 * t47 + (qJD(2) * t111 - t38) * t210) * t147, -t45 * t122 - t19 * t149 - t51 * t98 - t54 * t74 + (t122 * t203 - pkin(2) * t88 + t146 * t47 + (-qJD(2) * t112 + t39) * t210) * t147, -t111 * t88 - t112 * t87 + t38 * t98 + t39 * t97 + t44 * t74 + t45 * t73 + (-t146 * t18 + t148 * t19 + (t146 * t74 - t148 * t73) * qJD(3)) * t147, t18 * t111 + t19 * t112 - t38 * t44 - t39 * t45 - t51 * t54 + (-pkin(2) * t47 + (-t146 * t38 + t148 * t39) * qJD(3)) * t147, t57, -t16, -t243, t37, -t244, t36, -t26 * t98 - t27 * t97 - t41 * t74 + t93 * t87 + t94 * t88 + t229 * t73 + (-t148 * t9 + (qJD(3) * t74 + t14) * t146) * t147, t41 * t122 + t14 * t149 + t20 * t97 + t35 * t73 - t95 * t87 + (t148 * t5 + (qJD(4) * t73 - t205) * t146 + (qJD(2) * t94 - t26) * t210) * t147, -t9 * t149 + t20 * t98 + t35 * t74 - t95 * t88 + t229 * t122 + ((-t5 + t214) * t146 + (-qJD(2) * t93 + t27) * t210) * t147, t14 * t94 - t20 * t35 - t26 * t41 + t5 * t95 + t9 * t93 + t229 * t27 + (qJD(3) * t26 - qJD(4) * t20) * t225, t57, -t243, t16, t36, t244, t37, t13 * t97 - t6 * t98 + t63 * t88 - t68 * t87 + t230 * t74 - t228 * t73 + (t146 * t2 + t148 * t3) * t147, t3 * t149 - t69 * t88 + t7 * t98 - t231 * t74 - t228 * t122 + (-t1 * t146 + (qJD(2) * t68 - t13) * t210) * t147, -t2 * t149 + t69 * t87 - t7 * t97 + t231 * t73 + t230 * t122 + (-t1 * t148 + (-qJD(2) * t63 + t6) * t210) * t147, t1 * t69 + t228 * t13 + t2 * t63 + t230 * t6 + t231 * t7 + t3 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t162, t179, t38 * t52 - t39 * t53 + t47, 0, 0, 0, 0, 0, 0, t179, -t159 + t232, -t162, -t26 * t52 + t27 * t53 + t5, 0, 0, 0, 0, 0, 0, t179, -t162, t155, -t13 * t53 - t6 * t52 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, -t189, t247, -t27 * t122 + t20 * t74 + t14, 0, 0, 0, 0, 0, 0, t161, t247, t189, t7 * t74 + (qJD(5) + t13) * t122 + t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122 * t74 - t159, t136 + t239, -t73 ^ 2 - t119, -t6 * t122 - t7 * t73 + t3;];
tauc_reg = t17;
