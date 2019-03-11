% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:56
% EndTime: 2019-03-08 20:21:04
% DurationCPUTime: 2.85s
% Computational Cost: add. (3209->367), mult. (7170->484), div. (0->0), fcn. (4800->8), ass. (0->200)
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t106 = sin(qJ(4));
t104 = cos(pkin(6));
t109 = cos(qJ(4));
t215 = t104 * t109;
t170 = qJD(1) * t215;
t111 = -pkin(2) - pkin(8);
t110 = cos(qJ(2));
t103 = sin(pkin(6));
t207 = qJD(1) * t103;
t172 = t110 * t207;
t149 = qJD(3) - t172;
t72 = t111 * qJD(2) + t149;
t52 = t106 * t72 + t170;
t41 = qJD(4) * pkin(9) + t52;
t107 = sin(qJ(2));
t173 = t107 * t207;
t86 = t106 * pkin(4) - t109 * pkin(9) + qJ(3);
t61 = qJD(2) * t86 + t173;
t16 = t105 * t61 + t108 * t41;
t198 = qJD(5) * t108;
t199 = qJD(5) * t105;
t218 = t103 * t107;
t179 = qJD(2) * t218;
t201 = qJD(4) * t109;
t206 = qJD(1) * t106;
t33 = t72 * t201 + (-qJD(4) * t104 + t179) * t206;
t155 = pkin(4) * t109 + pkin(9) * t106;
t76 = qJD(4) * t155 + qJD(3);
t53 = (t76 + t172) * qJD(2);
t168 = t105 * t33 - t108 * t53 + t41 * t198 + t61 * t199;
t195 = t106 * qJD(2);
t97 = qJD(5) + t195;
t134 = t16 * t97 - t168;
t202 = qJD(4) * t106;
t193 = qJD(2) * qJD(4);
t214 = t105 * t106;
t204 = qJD(2) * t109;
t171 = t108 * t204;
t196 = t105 * qJD(4);
t81 = t171 + t196;
t57 = qJD(5) * t81 - t193 * t214;
t219 = t57 * t108;
t197 = qJD(5) * t109;
t176 = t105 * t197;
t194 = t108 * qJD(4);
t56 = qJD(2) * (t106 * t194 + t176) - qJD(5) * t194;
t220 = t56 * t105;
t229 = t108 * t81;
t79 = t105 * t204 - t194;
t267 = t109 * (qJD(5) * (t105 * t79 - t229) - t219 + t220) + (t105 * t81 + t108 * t79) * t202;
t266 = (t109 * (t79 + t194) - t97 * t214) * qJD(2) - t97 * t199;
t102 = t109 ^ 2;
t138 = qJD(2) * t102 - t106 * t97;
t174 = t108 * t197;
t224 = t109 * t79;
t265 = (t105 * t138 + t224) * qJD(4) + t106 * t57 + t97 * t174;
t226 = t109 * t56;
t227 = t108 * t97;
t264 = (t106 * (-t81 + t171) + t109 * t227) * qJD(4) - t226;
t135 = t81 * t97;
t136 = t79 * t97;
t263 = (t57 + t135) * t105 + (t56 + t136) * t108;
t12 = t97 * qJ(6) + t16;
t98 = t109 * t193;
t162 = pkin(5) * t98;
t2 = -t162 + t168;
t262 = -t12 * t97 + t2;
t200 = qJD(4) * t111;
t177 = t109 * t200;
t211 = t107 * t108;
t261 = t105 * t76 + t108 * t177 + t86 * t198 - (t105 * t110 + t106 * t211) * t207;
t259 = t57 - t135;
t51 = -t104 * t206 + t109 * t72;
t158 = qJD(2) * t173;
t85 = t109 * t158;
t34 = t52 * qJD(4) - t85;
t117 = -(t106 * t51 - t109 * t52) * qJD(4) + t33 * t106 - t34 * t109;
t212 = t106 * t111;
t237 = t105 * t86 + t108 * t212;
t252 = -qJD(5) * t237 + t173 * t214 + (-t172 + t76) * t108;
t251 = t81 ^ 2;
t250 = pkin(9) * t81;
t40 = -qJD(4) * pkin(4) - t51;
t17 = t79 * pkin(5) - t81 * qJ(6) + t40;
t247 = t17 * t81;
t217 = t103 * t110;
t70 = t104 * t106 + t109 * t217;
t246 = t34 * t70;
t5 = t57 * pkin(5) + t56 * qJ(6) - t81 * qJD(6) + t34;
t245 = t5 * t105;
t244 = t5 * t108;
t243 = t81 * t79;
t175 = t111 * t199;
t242 = qJ(6) * t201 + (qJD(6) - t175) * t106 + t261;
t169 = t105 * t111 - pkin(5);
t241 = t169 * t201 - t252;
t240 = -t106 * t175 + t261;
t239 = -t105 * t177 + t252;
t150 = pkin(5) * t105 - qJ(6) * t108;
t238 = t170 + (-qJD(2) * t150 + t72) * t106 - qJD(5) * t150 + t105 * qJD(6);
t83 = t155 * qJD(2);
t27 = t105 * t83 + t108 * t51;
t236 = qJD(2) * pkin(2);
t235 = t105 * t40;
t234 = t105 * t97;
t233 = t106 * t40;
t77 = (qJD(3) + t172) * qJD(2);
t232 = t107 * t77;
t231 = t108 * t40;
t228 = t108 * t86;
t205 = qJD(2) * qJ(3);
t84 = t173 + t205;
t223 = t110 * t84;
t222 = t34 * t105;
t221 = t34 * t108;
t113 = qJD(2) ^ 2;
t216 = t103 * t113;
t213 = t106 * t108;
t15 = -t105 * t41 + t108 * t61;
t210 = qJD(6) - t15;
t101 = t106 ^ 2;
t209 = t101 - t102;
t112 = qJD(4) ^ 2;
t208 = -t112 - t113;
t203 = qJD(2) * t110;
t192 = pkin(9) * t234;
t191 = pkin(9) * t227;
t189 = pkin(9) * t201;
t186 = t79 ^ 2 - t251;
t67 = t79 * t202;
t184 = t106 * t217;
t183 = t107 * t216;
t182 = t110 * t216;
t181 = t109 * t113 * t106;
t180 = t97 * t204;
t178 = t103 * t203;
t167 = qJD(5) * t79 - t56;
t166 = t97 + t195;
t164 = qJD(5) * t106 + qJD(2);
t163 = t81 * t173;
t161 = t109 * t173;
t160 = t109 * t179;
t159 = t106 * t179;
t157 = qJ(6) * t98;
t156 = t106 * t98;
t154 = -t84 + t173;
t152 = t167 * pkin(9);
t151 = t108 * pkin(5) + t105 * qJ(6);
t11 = -t97 * pkin(5) + t210;
t148 = t105 * t12 - t108 * t11;
t147 = t105 * t11 + t108 * t12;
t146 = t105 * t16 + t108 * t15;
t145 = t105 * t15 - t108 * t16;
t26 = -t105 * t51 + t108 * t83;
t139 = t77 * qJ(3) + t84 * qJD(3);
t137 = t164 * t97;
t133 = -t106 * t17 + t189;
t132 = -t189 + t233;
t131 = -t111 + t150;
t71 = -t184 + t215;
t130 = t103 * t211 - t71 * t105;
t47 = t105 * t218 + t71 * t108;
t44 = -qJD(4) * t70 + t159;
t8 = qJD(5) * t47 + t44 * t105 - t108 * t178;
t9 = qJD(5) * t130 + t105 * t178 + t44 * t108;
t128 = t130 * t56 - t47 * t57 - t9 * t79 + t8 * t81;
t127 = -t105 * t53 - t108 * t33 - t61 * t198 + t41 * t199;
t125 = t154 - t205;
t45 = -qJD(4) * t184 + t104 * t201 - t160;
t123 = t130 * t98 + t45 * t79 + t70 * t57 - t8 * t97;
t122 = t105 * t136 - t219;
t121 = qJD(2) * t149 - t111 * t112 + t77;
t120 = -t45 * t81 + t47 * t98 + t70 * t56 + t9 * t97;
t1 = t97 * qJD(6) - t127 + t157;
t119 = -qJD(5) * t148 + t1 * t108 + t2 * t105;
t118 = -qJD(5) * t146 + t105 * t168 - t108 * t127;
t116 = t79 * t174 + (t109 * t57 - t67) * t105;
t115 = t67 - t108 * t137 + (-t166 * t196 - t57) * t109;
t114 = -t57 * t213 - t194 * t224 + t164 * t229 + (qJD(2) * t79 + t167 * t106 + t81 * t201) * t105;
t87 = -pkin(4) - t151;
t78 = t149 - t236;
t65 = t131 * t109;
t64 = t166 * t201;
t62 = -t105 * t212 + t228;
t60 = t79 * t161;
t55 = t169 * t106 - t228;
t54 = t106 * qJ(6) + t237;
t48 = pkin(9) * t219;
t43 = t81 * pkin(5) + t79 * qJ(6);
t30 = t136 - t56;
t25 = (qJD(5) * t151 - qJD(6) * t108) * t109 - t131 * t202;
t22 = -pkin(5) * t204 - t26;
t21 = qJ(6) * t204 + t27;
t20 = t97 * t198 + (t97 * t213 + (-t81 + t196) * t109) * qJD(2);
t13 = t108 * t135 - t220;
t10 = -t81 * t176 + (-t81 * t202 - t226) * t108;
t7 = -t97 * t176 - t56 * t106 + (t108 * t138 + t109 * t81) * qJD(4);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183, -t182, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, t182 (t232 + (t223 + (t78 - t172) * t107) * qJD(2)) * t103, 0, 0, 0, 0, 0, 0, t106 * t182 + (-t45 + t160) * qJD(4), t109 * t182 + (-t44 - t159) * qJD(4) (-t106 * t44 + t109 * t45 + (-t106 * t70 - t109 * t71) * qJD(4)) * qJD(2), t33 * t71 + t246 + t52 * t44 - t51 * t45 + (t84 * t203 + t232) * t103, 0, 0, 0, 0, 0, 0, t123, -t120, t128, -t127 * t47 - t130 * t168 - t15 * t8 + t16 * t9 + t40 * t45 + t246, 0, 0, 0, 0, 0, 0, t123, t128, t120, t1 * t47 + t11 * t8 + t12 * t9 - t130 * t2 + t17 * t45 + t5 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3) (-t223 + (-t78 - t236) * t107) * t207 + t139, -0.2e1 * t156, 0.2e1 * t209 * t193, -t112 * t106, 0.2e1 * t156, -t112 * t109, 0, t106 * t121 - t125 * t201, t109 * t121 + t125 * t202 (t101 + t102) * t158 - t117 (-t223 + (-t106 * t52 - t109 * t51) * t107) * t207 + t117 * t111 + t139, t10, t267, t7, t116, -t265, t64, t60 + t239 * t97 + (-t168 + (t111 * t79 - t235) * qJD(4)) * t106 + (t40 * t198 + t222 - t111 * t57 + (qJD(2) * t62 + t15) * qJD(4)) * t109, -t240 * t97 + (t127 + (t111 * t81 - t231) * qJD(4)) * t106 + (t163 - t40 * t199 + t221 + t111 * t56 + (-qJD(2) * t237 - t16) * qJD(4)) * t109, t62 * t56 - t237 * t57 - t239 * t81 - t240 * t79 + t146 * t202 + (qJD(5) * t145 + t105 * t127 + t108 * t168) * t109, t200 * t233 - t127 * t237 - t168 * t62 + t240 * t16 + t239 * t15 + (-t111 * t34 + t173 * t40) * t109, t10, t7, -t267, t64, t265, t116, t25 * t79 + t65 * t57 + t60 - t241 * t97 + (-t17 * t196 - t2) * t106 + (t17 * t198 + t245 + (-qJD(2) * t55 - t11) * qJD(4)) * t109, -t54 * t57 - t55 * t56 + t241 * t81 - t242 * t79 + t148 * t202 + (-qJD(5) * t147 - t1 * t105 + t108 * t2) * t109, -t25 * t81 + t65 * t56 + t242 * t97 + (t17 * t194 + t1) * t106 + (-t163 + t17 * t199 - t244 + (qJD(2) * t54 + t12) * qJD(4)) * t109, t1 * t54 + t2 * t55 + t5 * t65 + (t25 + t161) * t17 + t242 * t12 + t241 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t154 * qJD(2), 0, 0, 0, 0, 0, 0, t208 * t106, t208 * t109, 0, -t84 * qJD(2) + t117, 0, 0, 0, 0, 0, 0, t115, t164 * t234 - t264, t114, -t146 * qJD(2) + (-qJD(4) * t145 - t34) * t109 + (qJD(4) * t40 + t118) * t106, 0, 0, 0, 0, 0, 0, t115, t114, -t105 * t137 + t264, -t148 * qJD(2) + (qJD(4) * t147 - t5) * t109 + (qJD(4) * t17 + t119) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, -t209 * t113, 0, -t181, 0, 0, -t84 * t204 + t85, -t154 * t195, 0, 0, t13, -t263, t20, t122, t266, -t180, -pkin(4) * t57 - t221 - t26 * t97 - t52 * t79 + (-t191 + t235) * qJD(5) + (t105 * t132 - t109 * t15) * qJD(2), pkin(4) * t56 + t222 + t27 * t97 - t52 * t81 + (t192 + t231) * qJD(5) + (t108 * t132 + t109 * t16) * qJD(2), t26 * t81 + t27 * t79 - t48 + (-t15 * t195 - t127 + (-t15 + t250) * qJD(5)) * t108 + (t152 - t134) * t105, -t34 * pkin(4) + pkin(9) * t118 - t15 * t26 - t16 * t27 - t40 * t52, t13, t20, t263, -t180, -t266, t122, -t244 + t22 * t97 + t87 * t57 - t238 * t79 + (t105 * t17 - t191) * qJD(5) + (-t105 * t133 + t109 * t11) * qJD(2), t21 * t79 - t22 * t81 - t48 + (t11 * t195 + t1 + (t11 + t250) * qJD(5)) * t108 + (t152 + t262) * t105, -t245 - t21 * t97 + t87 * t56 + t238 * t81 + (-t108 * t17 - t192) * qJD(5) + (t108 * t133 - t109 * t12) * qJD(2), pkin(9) * t119 - t11 * t22 - t12 * t21 - t17 * t238 + t5 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, -t186, t30, -t243, -t259, t98, -t40 * t81 + t134, t15 * t97 + t40 * t79 + t127, 0, 0, t243, t30, t186, t98, t259, -t243, -t43 * t79 + t134 + 0.2e1 * t162 - t247, pkin(5) * t56 - t57 * qJ(6) + (t12 - t16) * t81 + (t11 - t210) * t79, 0.2e1 * t157 - t17 * t79 + t43 * t81 + (0.2e1 * qJD(6) - t15) * t97 - t127, -t2 * pkin(5) + t1 * qJ(6) - t11 * t16 + t12 * t210 - t17 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98 + t243, t30, -t97 ^ 2 - t251, t247 + t262;];
tauc_reg  = t3;
