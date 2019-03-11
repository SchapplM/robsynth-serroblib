% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:09
% EndTime: 2019-03-08 18:58:17
% DurationCPUTime: 3.49s
% Computational Cost: add. (4931->385), mult. (13295->534), div. (0->0), fcn. (11347->12), ass. (0->199)
t128 = sin(qJ(5));
t278 = pkin(9) * t128;
t129 = sin(qJ(4));
t132 = cos(qJ(4));
t173 = pkin(4) * t129 - pkin(10) * t132;
t105 = t173 * qJD(4);
t124 = sin(pkin(6));
t122 = sin(pkin(12));
t130 = sin(qJ(3));
t133 = cos(qJ(3));
t125 = cos(pkin(12));
t126 = cos(pkin(7));
t221 = t125 * t126;
t155 = t122 * t133 + t130 * t221;
t144 = t155 * t124;
t127 = cos(pkin(6));
t112 = qJD(1) * t127 + qJD(2);
t123 = sin(pkin(7));
t224 = t123 * t130;
t196 = t112 * t224;
t70 = qJD(1) * t144 + t196;
t277 = pkin(9) * qJD(5) * t132 - t105 + t70;
t131 = cos(qJ(5));
t211 = qJD(4) * t128;
t214 = qJD(3) * t129;
t102 = t131 * t214 + t211;
t212 = qJD(3) * t132;
t113 = -qJD(5) + t212;
t226 = t102 * t113;
t205 = t131 * qJD(4);
t100 = t128 * t214 - t205;
t229 = t100 * t113;
t203 = qJD(4) * qJD(5);
t208 = qJD(5) * t128;
t186 = t129 * t208;
t273 = -t132 * t205 + t186;
t84 = t273 * qJD(3) - t131 * t203;
t207 = qJD(5) * t131;
t209 = qJD(4) * t132;
t148 = t128 * t209 + t129 * t207;
t85 = qJD(3) * t148 + t128 * t203;
t276 = (t84 - t229) * t131 + (t85 - t226) * t128;
t215 = qJD(1) * t124;
t191 = t125 * t215;
t86 = t112 * t126 - t123 * t191;
t238 = t129 * t86;
t66 = qJD(3) * pkin(9) + t70;
t41 = t132 * t66 + t238;
t39 = qJD(4) * pkin(10) + t41;
t110 = -pkin(4) * t132 - pkin(10) * t129 - pkin(3);
t225 = t122 * t130;
t69 = t133 * (t112 * t123 + t126 * t191) - t215 * t225;
t53 = qJD(3) * t110 - t69;
t12 = t128 * t53 + t131 * t39;
t266 = -t129 * t66 + t132 * t86;
t223 = t123 * t133;
t264 = (t133 * t221 - t225) * t124;
t270 = qJD(1) * t264 + t112 * t223;
t61 = t270 * qJD(3);
t15 = qJD(4) * t266 + t132 * t61;
t49 = (t105 + t70) * qJD(3);
t181 = t128 * t15 - t131 * t49 + t39 * t207 + t53 * t208;
t157 = -t113 * t12 - t181;
t272 = t127 * t223 + t264;
t230 = t85 * t131;
t231 = t84 * t128;
t271 = t129 * (qJD(5) * (t100 * t128 - t102 * t131) - t230 + t231) - (t100 * t131 + t102 * t128) * t209;
t220 = t128 * t132;
t269 = qJD(3) * (t129 * (t100 + t205) - t113 * t220) + t113 * t208;
t120 = t129 ^ 2;
t163 = qJD(3) * t120 - t113 * t132;
t185 = t113 * t207;
t268 = qJD(4) * (t100 * t129 + t128 * t163) - t129 * t185 - t132 * t85;
t10 = -qJ(6) * t113 + t12;
t204 = qJD(3) * qJD(4);
t117 = t129 * t204;
t178 = pkin(5) * t117;
t2 = -t178 + t181;
t267 = t10 * t113 + t2;
t218 = t131 * t132;
t247 = t129 * t205 * pkin(9) - t110 * t207 + t277 * t128 + t69 * t218;
t265 = t85 + t226;
t210 = qJD(4) * t129;
t249 = -t110 * t208 - t277 * t131 + t210 * t278 + t69 * t220;
t259 = t102 ^ 2;
t258 = pkin(10) * t102;
t257 = pkin(10) * t113;
t16 = t41 * qJD(4) + t129 * t61;
t8 = t85 * pkin(5) + t84 * qJ(6) - t102 * qJD(6) + t16;
t256 = t128 * t8;
t255 = t131 * t8;
t75 = t127 * t224 + t144;
t93 = -t123 * t124 * t125 + t126 * t127;
t50 = t129 * t75 - t93 * t132;
t254 = t16 * t50;
t94 = -t132 * t126 + t129 * t224;
t253 = t16 * t94;
t62 = t70 * qJD(3);
t252 = t62 * t272;
t169 = pkin(5) * t128 - qJ(6) * t131;
t251 = t238 + (t169 * qJD(3) + t66) * t132 - t169 * qJD(5) + t128 * qJD(6);
t250 = -pkin(5) * t210 - t249;
t248 = -qJ(6) * t210 + t132 * qJD(6) + t247;
t104 = t173 * qJD(3);
t27 = t128 * t104 + t131 * t266;
t246 = qJD(3) * pkin(3);
t244 = t100 * t69;
t38 = -qJD(4) * pkin(4) - t266;
t19 = pkin(5) * t100 - qJ(6) * t102 + t38;
t243 = t102 * t19;
t242 = t102 * t69;
t240 = t128 * t38;
t239 = t129 * t69;
t237 = t131 * t38;
t235 = t133 * t62;
t234 = t16 * t128;
t233 = t16 * t129;
t232 = t16 * t131;
t89 = pkin(9) * t218 + t128 * t110;
t227 = t102 * t100;
t135 = qJD(3) ^ 2;
t222 = t123 * t135;
t219 = t131 * t110;
t11 = -t128 * t39 + t131 * t53;
t217 = qJD(6) - t11;
t121 = t132 ^ 2;
t216 = t120 - t121;
t213 = qJD(3) * t130;
t202 = t128 * t257;
t201 = t131 * t257;
t200 = pkin(10) * t210;
t199 = pkin(10) * t205;
t197 = t128 * t223;
t194 = t130 * t222;
t192 = t129 * t135 * t132;
t190 = t123 * t213;
t189 = qJD(3) * t223;
t184 = t113 * t214;
t183 = t100 ^ 2 - t259;
t65 = -t69 - t246;
t180 = -qJD(3) * t65 - t61;
t177 = t132 * t189;
t176 = t129 * t189;
t175 = qJ(6) * t117;
t174 = t132 * t117;
t9 = pkin(5) * t113 + t217;
t172 = -t10 * t128 + t131 * t9;
t171 = (qJD(5) * t100 - t84) * pkin(10);
t170 = pkin(5) * t131 + qJ(6) * t128;
t168 = -t11 * t131 - t12 * t128;
t51 = t129 * t93 + t132 * t75;
t25 = -t128 * t272 + t131 * t51;
t24 = t128 * t51 + t131 * t272;
t26 = t104 * t131 - t128 * t266;
t160 = pkin(9) + t169;
t134 = qJD(4) ^ 2;
t159 = pkin(9) * t134;
t158 = qJD(4) * (t65 + t69 - t246);
t95 = t126 * t129 + t132 * t224;
t79 = t128 * t95 + t131 * t223;
t154 = -t128 * t49 - t131 * t15 - t53 * t207 + t39 * t208;
t71 = t272 * qJD(3);
t21 = -qJD(4) * t50 + t71 * t132;
t72 = t75 * qJD(3);
t6 = qJD(5) * t25 + t21 * t128 - t72 * t131;
t7 = -qJD(5) * t24 + t72 * t128 + t21 * t131;
t153 = -t7 * t100 + t102 * t6 - t24 * t84 - t25 * t85;
t77 = -qJD(4) * t94 + t177;
t36 = -qJD(5) * t79 + t128 * t190 + t131 * t77;
t37 = -qJD(5) * t197 + t128 * t77 - t131 * t190 + t95 * t207;
t80 = t131 * t95 - t197;
t151 = -t36 * t100 + t102 * t37 - t79 * t84 - t80 * t85;
t20 = qJD(4) * t51 + t71 * t129;
t143 = t20 * t100 + t113 * t6 - t24 * t117 + t50 * t85;
t78 = qJD(4) * t95 + t176;
t142 = t78 * t100 + t113 * t37 - t79 * t117 + t94 * t85;
t141 = -t128 * t229 - t230;
t139 = -t102 * t20 - t113 * t7 + t25 * t117 + t50 * t84;
t138 = -t102 * t78 - t113 * t36 + t80 * t117 + t84 * t94;
t137 = t233 + t15 * t132 + (-t129 * t41 - t132 * t266) * qJD(4);
t136 = t85 * t128 * t129 + t100 * t148;
t107 = -pkin(4) - t170;
t91 = t160 * t129;
t88 = -pkin(9) * t220 + t219;
t87 = (-t113 - t212) * t210;
t83 = -t219 + (pkin(5) + t278) * t132;
t82 = -qJ(6) * t132 + t89;
t81 = pkin(10) * t230;
t76 = pkin(5) * t102 + qJ(6) * t100;
t60 = -t84 - t229;
t56 = (t170 * qJD(5) - qJD(6) * t131) * t129 + t160 * t209;
t55 = -t185 + (t113 * t218 + (-t102 + t211) * t129) * qJD(3);
t43 = -t131 * t226 - t231;
t42 = -t84 * t131 * t129 - t273 * t102;
t32 = t113 * t186 + t84 * t132 + (t102 * t129 + t131 * t163) * qJD(4);
t23 = -pkin(5) * t214 - t26;
t22 = qJ(6) * t214 + t27;
t1 = -qJD(6) * t113 - t154 + t175;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 * qJD(3), -t71 * qJD(3), 0, t61 * t75 - t69 * t72 + t70 * t71 - t252, 0, 0, 0, 0, 0, 0, -t20 * qJD(4) + (-t132 * t72 - t210 * t272) * qJD(3), -t21 * qJD(4) + (t129 * t72 - t209 * t272) * qJD(3) (t129 * t20 + t132 * t21 + (-t129 * t51 + t132 * t50) * qJD(4)) * qJD(3), t15 * t51 - t20 * t266 + t21 * t41 + t65 * t72 - t252 + t254, 0, 0, 0, 0, 0, 0, t143, -t139, t153, -t11 * t6 + t12 * t7 - t154 * t25 + t181 * t24 + t20 * t38 + t254, 0, 0, 0, 0, 0, 0, t143, t153, t139, t1 * t25 + t10 * t7 + t19 * t20 + t2 * t24 + t50 * t8 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t194, -t133 * t222, 0 (t130 * t61 - t235 + (-t130 * t69 + t133 * t70) * qJD(3)) * t123, 0, 0, 0, 0, 0, 0, -t132 * t194 + (-t78 - t176) * qJD(4), t129 * t194 + (-t77 - t177) * qJD(4) (t129 * t78 + t132 * t77 + (-t129 * t95 + t132 * t94) * qJD(4)) * qJD(3), t15 * t95 + t253 - t266 * t78 + t41 * t77 + (t213 * t65 - t235) * t123, 0, 0, 0, 0, 0, 0, t142, -t138, t151, -t11 * t37 + t12 * t36 - t154 * t80 + t181 * t79 + t38 * t78 + t253, 0, 0, 0, 0, 0, 0, t142, t151, t138, t1 * t80 + t10 * t36 + t19 * t78 + t2 * t79 + t37 * t9 + t8 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t155 * t215 - t196 + t70) * qJD(3) (-t270 + t69) * qJD(3), 0, 0, 0.2e1 * t174, -0.2e1 * t216 * t204, t134 * t132, -0.2e1 * t174, -t134 * t129, 0, t129 * t158 - t132 * t159, t129 * t159 + t132 * t158 (-t120 - t121) * t69 * qJD(3) + t137, -t62 * pkin(3) - t65 * t70 + (t129 * t266 - t132 * t41) * t69 + t137 * pkin(9), t42, t271, t32, t136, -t268, t87, -t249 * t113 + (t181 + (pkin(9) * t100 + t240) * qJD(4)) * t132 + (t38 * t207 + pkin(9) * t85 - t244 + t234 + (qJD(3) * t88 + t11) * qJD(4)) * t129, -t247 * t113 + (-t154 + (pkin(9) * t102 + t237) * qJD(4)) * t132 + (-t38 * t208 - pkin(9) * t84 - t242 + t232 + (-qJD(3) * t89 - t12) * qJD(4)) * t129, t84 * t88 - t85 * t89 - t249 * t102 + t247 * t100 + t168 * t209 + (t128 * t154 + t131 * t181 + (t11 * t128 - t12 * t131) * qJD(5)) * t129, -t38 * t239 - t154 * t89 - t181 * t88 - t247 * t12 + t249 * t11 + (t209 * t38 + t233) * pkin(9), t42, t32, -t271, t87, t268, t136, t100 * t56 + t85 * t91 + (t19 * t211 + t2) * t132 + t250 * t113 + (t19 * t207 - t244 + t256 + (-qJD(3) * t83 - t9) * qJD(4)) * t129, -t82 * t85 - t83 * t84 + t250 * t102 + t248 * t100 + t172 * t209 + (-t1 * t128 + t131 * t2 + (-t10 * t131 - t128 * t9) * qJD(5)) * t129, -t102 * t56 + t84 * t91 + (-t19 * t205 - t1) * t132 + t248 * t113 + (t19 * t208 + t242 - t255 + (qJD(3) * t82 + t10) * qJD(4)) * t129, t1 * t82 + t2 * t83 + t8 * t91 + t250 * t9 + (t56 - t239) * t19 - t248 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, t216 * t135, 0, t192, 0, 0, t180 * t129, t180 * t132, 0, 0, t43, -t276, t55, t141, t269, t184, -pkin(4) * t85 - t41 * t100 + t26 * t113 - t232 + (t201 + t240) * qJD(5) + (-t11 * t129 + (-t132 * t38 - t200) * t128) * qJD(3), pkin(4) * t84 - t41 * t102 - t27 * t113 + t234 + (-t202 + t237) * qJD(5) + (-t38 * t218 + (t12 - t199) * t129) * qJD(3), t100 * t27 + t102 * t26 - t81 + (t11 * t212 - t154 + (-t11 + t258) * qJD(5)) * t131 + (t171 - t157) * t128, -t16 * pkin(4) - t11 * t26 - t12 * t27 - t38 * t41 + (qJD(5) * t168 + t128 * t181 - t131 * t154) * pkin(10), t43, t55, t276, t184, -t269, t141, t107 * t85 - t113 * t23 - t255 - t251 * t100 + (t128 * t19 + t201) * qJD(5) + (t129 * t9 + (-t132 * t19 - t200) * t128) * qJD(3), t100 * t22 - t102 * t23 - t81 + (-t9 * t212 + t1 + (t9 + t258) * qJD(5)) * t131 + (t171 + t267) * t128, t107 * t84 + t113 * t22 - t256 + t251 * t102 + (-t131 * t19 + t202) * qJD(5) + (t19 * t218 + (-t10 + t199) * t129) * qJD(3), -t10 * t22 + t8 * t107 - t9 * t23 - t251 * t19 + (qJD(5) * t172 + t1 * t131 + t2 * t128) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, -t183, t60, -t227, -t265, t117, -t102 * t38 + t157, t100 * t38 - t11 * t113 + t154, 0, 0, t227, t60, t183, t117, t265, -t227, -t100 * t76 + t157 + 0.2e1 * t178 - t243, pkin(5) * t84 - t85 * qJ(6) + (t10 - t12) * t102 + (t9 - t217) * t100, 0.2e1 * t175 - t100 * t19 + t102 * t76 + (-0.2e1 * qJD(6) + t11) * t113 - t154, -t2 * pkin(5) + t1 * qJ(6) + t10 * t217 - t9 * t12 - t19 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117 + t227, t60, -t113 ^ 2 - t259, t243 + t267;];
tauc_reg  = t3;
