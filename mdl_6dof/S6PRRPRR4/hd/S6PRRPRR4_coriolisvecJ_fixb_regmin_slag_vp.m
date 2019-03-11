% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:28
% EndTime: 2019-03-08 22:14:35
% DurationCPUTime: 2.32s
% Computational Cost: add. (2301->306), mult. (5570->420), div. (0->0), fcn. (4080->10), ass. (0->186)
t127 = sin(qJ(5));
t128 = sin(qJ(3));
t131 = cos(qJ(5));
t132 = cos(qJ(3));
t78 = t127 * t128 + t131 * t132;
t242 = t78 * qJD(2);
t249 = qJD(6) + t242;
t118 = qJD(3) - qJD(5);
t126 = sin(qJ(6));
t130 = cos(qJ(6));
t199 = qJD(2) * t132;
t179 = t127 * t199;
t201 = qJD(2) * t128;
t74 = t131 * t201 - t179;
t55 = t130 * t118 + t126 * t74;
t254 = t249 * t55;
t154 = t118 * t126 - t130 * t74;
t253 = t154 * t249;
t124 = cos(pkin(6));
t202 = qJD(1) * t124;
t102 = t132 * t202;
t129 = sin(qJ(2));
t123 = sin(pkin(6));
t203 = qJD(1) * t123;
t182 = t129 * t203;
t88 = qJD(2) * pkin(8) + t182;
t75 = t128 * t88;
t60 = -t75 + t102;
t250 = qJD(4) - t60;
t252 = qJD(6) - t249;
t251 = -pkin(9) * t201 + t250;
t195 = qJD(5) * t131;
t197 = qJD(3) * t132;
t144 = t127 * t197 + t128 * t195;
t192 = qJD(2) * qJD(3);
t177 = t128 * t192;
t40 = qJD(2) * t144 - qJD(5) * t179 - t131 * t177;
t216 = t126 * t40;
t243 = t130 * t249;
t248 = t154 * t74 + t243 * t249 + t216;
t193 = qJD(6) * t130;
t194 = qJD(6) * t126;
t147 = t78 * qJD(5);
t176 = t132 * t192;
t221 = t127 * t177 + t131 * t176;
t39 = -qJD(2) * t147 + t221;
t14 = -t118 * t193 + t130 * t39 - t74 * t194;
t213 = t14 * t126;
t247 = t154 * t243 - t213;
t15 = -t154 * qJD(6) + t126 * t39;
t246 = (-t15 + t253) * t126 - (-t14 + t254) * t130;
t134 = -pkin(3) - pkin(4);
t183 = t134 * qJD(3);
t38 = t183 + t251;
t120 = qJD(3) * qJ(4);
t61 = t128 * t202 + t132 * t88;
t49 = -pkin(9) * t199 + t61;
t42 = t120 + t49;
t12 = -t127 * t42 + t131 * t38;
t10 = pkin(5) * t118 - t12;
t245 = t249 * t10;
t54 = t120 + t61;
t51 = -qJD(3) * pkin(3) + t250;
t244 = t118 * t74 + t40;
t164 = t128 * t183;
t114 = t128 * qJD(4);
t205 = qJ(4) * t197 + t114;
t159 = t164 + t205 + t182;
t196 = qJD(5) * t127;
t119 = qJD(3) * qJD(4);
t133 = cos(qJ(2));
t101 = t133 * t203;
t161 = qJD(2) * t101;
t222 = qJD(3) * t102 + t132 * t161;
t191 = t119 + t222;
t198 = qJD(3) * t128;
t25 = (pkin(9) * qJD(2) - t88) * t198 + t191;
t178 = qJD(1) * t198;
t35 = t124 * t178 + t128 * t161 + t88 * t197;
t27 = -pkin(9) * t176 + t35;
t2 = t127 * t25 - t131 * t27 + t42 * t195 + t38 * t196;
t41 = pkin(5) * t74 + pkin(10) * t242;
t110 = qJ(4) * t199;
t66 = t134 * t201 + t110;
t153 = qJ(4) * t131 + t127 * t134;
t82 = -pkin(10) + t153;
t241 = (qJD(6) * t82 - t41 + t66) * t249 - t2;
t240 = (pkin(10) * qJD(6) + t41) * t249 + t2;
t1 = t127 * t27 + t131 * t25 + t38 * t195 - t42 * t196;
t151 = t127 * t132 - t128 * t131;
t125 = qJD(2) * pkin(2);
t89 = -t101 - t125;
t62 = -pkin(3) * t199 - qJ(4) * t201 + t89;
t50 = pkin(4) * t199 - t62;
t18 = pkin(5) * t242 - pkin(10) * t74 + t50;
t234 = pkin(8) - pkin(9);
t91 = t234 * t128;
t92 = t234 * t132;
t58 = t127 * t92 - t131 * t91;
t83 = t234 * t198;
t84 = qJD(3) * t92;
t225 = -t58 * qJD(5) - t78 * t101 + t127 * t84 - t131 * t83;
t90 = -t132 * pkin(3) - t128 * qJ(4) - pkin(2);
t77 = t132 * pkin(4) - t90;
t28 = pkin(5) * t78 + pkin(10) * t151 + t77;
t44 = t78 * qJD(3) - t147;
t59 = t127 * t91 + t131 * t92;
t237 = -(qJD(6) * t18 + t1) * t78 + t10 * t44 - t2 * t151 + (-qJD(6) * t28 - t225) * t249 - t59 * t40;
t214 = t130 * t40;
t236 = t126 * t249 ^ 2 - t55 * t74 - t214;
t235 = -(t60 + t75) * qJD(3) + t222;
t13 = t127 * t38 + t131 * t42;
t11 = -pkin(10) * t118 + t13;
t155 = t11 * t126 - t130 * t18;
t233 = t155 * t74;
t4 = t11 * t130 + t126 * t18;
t232 = t4 * t74;
t231 = t10 * t151;
t230 = t28 * t40;
t229 = t40 * t151;
t228 = t249 * t74;
t227 = t74 * t242;
t152 = -qJ(4) * t127 + t131 * t134;
t226 = -t152 * qJD(5) + t127 * t49 - t251 * t131;
t224 = t59 * qJD(5) - t151 * t101 - t127 * t83 - t131 * t84;
t223 = t153 * qJD(5) + t251 * t127 + t131 * t49;
t219 = t118 * t242;
t212 = t123 * t129;
t211 = t123 * t133;
t136 = qJD(2) ^ 2;
t210 = t123 * t136;
t135 = qJD(3) ^ 2;
t209 = t135 * t128;
t208 = t135 * t132;
t206 = -qJ(4) * t176 - qJD(2) * t114;
t121 = t128 ^ 2;
t122 = t132 ^ 2;
t204 = t121 - t122;
t200 = qJD(2) * t129;
t190 = pkin(3) * t198;
t189 = t242 ^ 2 - t74 ^ 2;
t188 = t249 * t201;
t187 = t151 * t194;
t186 = t249 * t193;
t185 = t129 * t210;
t184 = t128 * t136 * t132;
t181 = t123 * t200;
t180 = qJD(2) * t211;
t175 = t89 - t125;
t171 = qJD(2) * t90 + t62;
t166 = t118 * t249;
t165 = t118 ^ 2;
t163 = t128 * t180;
t162 = t132 * t180;
t45 = -t131 * t198 - t132 * t196 + t144;
t160 = pkin(5) * t45 - pkin(10) * t44 + t159;
t158 = t249 * t44 - t229;
t70 = -t124 * t132 + t128 * t212;
t71 = t124 * t128 + t132 * t212;
t32 = t127 * t71 - t131 * t70;
t33 = t127 * t70 + t131 * t71;
t150 = qJD(3) * t61 - t35;
t149 = -t126 * t33 + t130 * t211;
t148 = t126 * t211 + t130 * t33;
t145 = t50 * t74 + t2;
t143 = -pkin(10) * t40 + t12 * t249 + t245;
t142 = -t242 * t50 + t1;
t43 = (t182 + t190) * qJD(2) + t206;
t67 = t190 - t205;
t141 = -pkin(8) * t135 - t43 + (-t67 + t182) * qJD(2);
t140 = t226 * t249 - t82 * t40 - t245;
t29 = -t88 * t198 + t191;
t139 = t128 * t35 + t132 * t29 + (-t128 * t54 + t132 * t51) * qJD(3);
t47 = qJD(3) * t71 + t163;
t138 = -t132 * t185 + (-t47 - t163) * qJD(3);
t36 = (t164 - t182) * qJD(2) - t206;
t85 = t178 * t211;
t81 = pkin(5) - t152;
t80 = pkin(3) * t201 - t110;
t46 = -qJD(3) * t70 + t162;
t21 = -t128 * t185 + (t46 + t162) * qJD(3);
t8 = -t32 * qJD(5) + t127 * t47 + t131 * t46;
t7 = t33 * qJD(5) + t127 * t46 - t131 * t47;
t6 = pkin(5) * t40 - pkin(10) * t39 + t36;
t5 = t130 * t6;
t3 = [0, 0, -t185, -t133 * t210, 0, 0, 0, 0, 0, t138, -t21, t138 (t128 * t47 + t132 * t46 + (-t128 * t71 + t132 * t70) * qJD(3)) * qJD(2), t21, t29 * t71 + t35 * t70 + t46 * t54 + t47 * t51 + (-t133 * t43 + t62 * t200) * t123, 0, 0, 0, 0, 0, t118 * t7 + (t133 * t40 - t200 * t242) * t123, t118 * t8 + (t133 * t39 - t74 * t200) * t123, 0, 0, 0, 0, 0 (-qJD(6) * t148 - t126 * t8 - t130 * t181) * t249 + t149 * t40 + t7 * t55 + t32 * t15 -(qJD(6) * t149 - t126 * t181 + t130 * t8) * t249 - t148 * t40 - t7 * t154 + t32 * t14; 0, 0, 0, 0, 0.2e1 * t128 * t176, -0.2e1 * t204 * t192, t208, -t209, 0, -pkin(8) * t208 + t175 * t198 + t85, pkin(8) * t209 + (t175 + t101) * t197, t141 * t132 + t171 * t198 + t85 (-t121 - t122) * t161 + t139 (-t171 - t101) * t197 + t141 * t128, t43 * t90 + t62 * t67 + (-t129 * t62 + (-t128 * t51 - t132 * t54) * t133) * t203 + t139 * pkin(8), -t151 * t39 + t44 * t74, -t242 * t44 - t39 * t78 - t45 * t74 + t229, -t44 * t118, t45 * t118, 0, t224 * t118 + t159 * t242 + t36 * t78 + t40 * t77 + t45 * t50, t225 * t118 - t151 * t36 + t159 * t74 + t39 * t77 + t44 * t50, -t154 * t187 + (-t14 * t151 - t154 * t44) * t130 (t126 * t154 - t130 * t55) * t44 - (-t213 - t130 * t15 + (t126 * t55 + t130 * t154) * qJD(6)) * t151, t158 * t130 + t14 * t78 - t154 * t45 + t187 * t249, -t158 * t126 - t15 * t78 + t151 * t186 - t45 * t55, t249 * t45 + t40 * t78, t58 * t15 - t155 * t45 + t5 * t78 + t224 * t55 + (t230 + t160 * t249 + (-t11 * t78 - t249 * t59 - t231) * qJD(6)) * t130 + t237 * t126, t58 * t14 - t4 * t45 - t224 * t154 + (-t230 - (-qJD(6) * t11 + t6) * t78 + qJD(6) * t231 + (qJD(6) * t59 - t160) * t249) * t126 + t237 * t130; 0, 0, 0, 0, -t184, t204 * t136, 0, 0, 0, -t89 * t201 + t150, -t89 * t199 - t235 (-t128 * t62 + t132 * t80) * qJD(2) + t150, 0, 0.2e1 * t119 + (t128 * t80 + t132 * t62) * qJD(2) + t235, -pkin(3) * t35 + qJ(4) * t29 + t250 * t54 - t51 * t61 - t62 * t80, -t227, t189, qJD(5) * t242 + t219 - t221, t244, 0, t223 * t118 - t242 * t66 + t145, -t226 * t118 - t66 * t74 + t142, t247, -t246, -t248, t236, t228, t140 * t126 - t241 * t130 + t81 * t15 + t223 * t55 - t233, t241 * t126 + t140 * t130 + t81 * t14 - t154 * t223 - t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, 0, -t121 * t136 - t135, -qJD(3) * t54 + t62 * t201 + t35, 0, 0, 0, 0, 0, -t127 * t165 - t201 * t242, -t131 * t165 - t74 * t201, 0, 0, 0, 0, 0, -t130 * t188 + (t126 * t166 - t15) * t131 + (-t118 * t55 - t186 - t216) * t127, t126 * t188 + (t130 * t166 - t14) * t131 + (t118 * t154 + t194 * t249 - t214) * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, -t189, t39 - t219, -t244, 0, -t118 * t13 - t145, -t118 * t12 - t142, -t247, t246, t248, -t236, -t228, -pkin(5) * t15 + t143 * t126 - t13 * t55 - t240 * t130 + t233, -pkin(5) * t14 + t240 * t126 + t13 * t154 + t143 * t130 + t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154 * t55, t154 ^ 2 - t55 ^ 2, t14 + t254, -t15 - t253, t40, -t1 * t126 + t10 * t154 - t252 * t4 + t5, -t1 * t130 + t10 * t55 - t126 * t6 + t252 * t155;];
tauc_reg  = t3;
