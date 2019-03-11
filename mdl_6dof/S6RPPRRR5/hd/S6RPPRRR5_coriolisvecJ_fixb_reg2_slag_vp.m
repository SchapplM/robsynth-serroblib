% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:00
% EndTime: 2019-03-09 02:29:08
% DurationCPUTime: 2.69s
% Computational Cost: add. (4664->322), mult. (9346->429), div. (0->0), fcn. (5862->6), ass. (0->195)
t111 = sin(qJ(4));
t113 = cos(qJ(5));
t114 = cos(qJ(4));
t215 = sin(qJ(5));
t72 = t113 * t111 + t215 * t114;
t67 = t72 * qJD(1);
t228 = qJD(6) + t67;
t110 = sin(qJ(6));
t112 = cos(qJ(6));
t175 = qJD(4) + qJD(5);
t151 = t112 * t175;
t184 = qJD(1) * t114;
t166 = t113 * t184;
t167 = t215 * t111;
t69 = -qJD(1) * t167 + t166;
t56 = t110 * t69 - t151;
t224 = t228 * t56;
t180 = qJD(6) * t110;
t176 = qJD(1) * qJD(4);
t162 = t111 * t176;
t181 = qJD(5) * t113;
t185 = qJD(1) * t111;
t163 = t215 * qJD(5);
t231 = qJD(4) * t215 + t163;
t227 = t231 * qJD(1);
t45 = t113 * t162 + t227 * t114 + t181 * t185;
t29 = -qJD(6) * t151 + t112 * t45 + t69 * t180;
t232 = -t29 + t224;
t100 = qJD(1) * qJ(2) + qJD(3);
t84 = -pkin(7) * qJD(1) + t100;
t160 = pkin(8) * qJD(1) - t84;
t64 = t160 * t111;
t194 = t113 * t64;
t65 = -pkin(8) * t184 + t114 * t84;
t62 = qJD(4) * pkin(4) + t65;
t40 = t215 * t62 - t194;
t36 = t175 * pkin(9) + t40;
t109 = pkin(1) + qJ(3);
t220 = qJD(1) * t109;
t85 = -qJD(2) + t220;
t71 = pkin(4) * t185 + t85;
t38 = t67 * pkin(5) - t69 * pkin(9) + t71;
t17 = t110 * t38 + t112 * t36;
t183 = qJD(4) * t111;
t177 = qJD(1) * qJD(2);
t94 = t114 * t177;
t126 = t160 * t183 + t94;
t119 = t215 * t126 + t64 * t163;
t178 = t111 * qJD(2);
t182 = qJD(4) * t114;
t55 = t84 * t182 + (-pkin(8) * t182 + t178) * qJD(1);
t158 = qJD(5) * t62 + t55;
t221 = t158 * t113;
t14 = t119 + t221;
t190 = t113 * t114;
t134 = t175 * t190;
t203 = t227 * t111;
t46 = qJD(1) * t134 - t203;
t104 = qJD(3) * qJD(1);
t161 = t114 * t176;
t75 = pkin(4) * t161 + t104;
t23 = t46 * pkin(5) + t45 * pkin(9) + t75;
t3 = -qJD(6) * t17 - t110 * t14 + t112 * t23;
t230 = -t17 * t228 - t3;
t130 = t110 * t36 - t112 * t38;
t2 = -qJD(6) * t130 + t110 * t23 + t112 * t14;
t229 = -t130 * t228 - t2;
t157 = t110 * t228;
t197 = t112 * t46;
t226 = t157 * t228 - t197;
t52 = -t231 * t111 + t134;
t223 = t52 * t175;
t222 = t67 * t175;
t106 = t111 ^ 2;
t107 = t114 ^ 2;
t187 = t106 + t107;
t219 = t187 * qJD(2);
t218 = t67 ^ 2;
t217 = t69 ^ 2;
t216 = 0.2e1 * t104;
t146 = -t113 * t126 + t215 * t55;
t15 = qJD(5) * t40 + t146;
t108 = -pkin(7) + qJ(2);
t205 = pkin(8) - t108;
t76 = t205 * t111;
t77 = t205 * t114;
t53 = t113 * t77 - t215 * t76;
t214 = t15 * t53;
t73 = -t167 + t190;
t213 = t15 * t73;
t1 = t2 * t112;
t212 = t46 * t72;
t211 = t56 * t69;
t58 = t110 * t175 + t112 * t69;
t210 = t58 * t56;
t209 = t58 * t69;
t208 = t228 * t69;
t207 = t69 * t67;
t206 = t73 * t46;
t41 = t215 * t65 - t194;
t204 = t40 - t41;
t202 = pkin(4) * qJD(5);
t201 = t110 * t45;
t200 = t110 * t46;
t199 = t110 * t56;
t198 = t110 * t67;
t196 = t112 * t58;
t195 = t112 * t67;
t193 = t29 * t110;
t30 = qJD(6) * t58 - t201;
t192 = t30 * t112;
t116 = qJD(1) ^ 2;
t191 = t106 * t116;
t189 = qJD(2) - t85;
t87 = pkin(4) * t182 + qJD(3);
t115 = qJD(4) ^ 2;
t186 = -t115 - t116;
t179 = qJD(6) * t112;
t174 = t130 * t195 - t17 * t198 + t1;
t91 = t111 * pkin(4) + t109;
t173 = pkin(4) * t181;
t172 = pkin(4) * t184;
t102 = 0.2e1 * t177;
t171 = t215 * t64;
t170 = t73 * t180;
t169 = t73 * t179;
t168 = t114 * t116 * t111;
t39 = t113 * t62 + t171;
t35 = -t175 * pkin(5) - t39;
t31 = t35 * t180;
t32 = t35 * t179;
t165 = t130 * t69 + t31;
t159 = t228 * t58;
t156 = t112 * t228;
t155 = t85 + t220;
t154 = qJD(1) * t187;
t153 = t69 * t175;
t152 = qJD(6) * t72 + qJD(1);
t150 = t15 * t110 + t17 * t69 + t32;
t149 = t111 * t161;
t48 = t69 * pkin(5) + t67 * pkin(9);
t145 = t130 * t52 - t3 * t72;
t144 = t17 * t52 + t2 * t72;
t51 = t175 * t72;
t143 = -t35 * t51 + t213;
t142 = -t73 * t29 - t51 * t58;
t141 = -t29 * t72 + t58 * t52;
t140 = t73 * t30 - t51 * t56;
t139 = -t30 * t72 - t56 * t52;
t138 = -t73 * t45 - t69 * t51;
t137 = -t228 * t51 + t206;
t136 = t52 * t67 + t212;
t135 = qJD(2) + t155;
t132 = t110 * t17 - t112 * t130;
t131 = -t110 * t130 - t112 * t17;
t47 = t72 * pkin(5) - t73 * pkin(9) + t91;
t54 = -t113 * t76 - t215 * t77;
t26 = -t110 * t54 + t112 * t47;
t27 = t110 * t47 + t112 * t54;
t129 = pkin(4) * t163 - t41;
t128 = -t108 * t115 + t216;
t127 = -t71 * t69 - t146;
t125 = t156 * t228 + t200;
t124 = t114 * qJD(2) + t205 * t183;
t123 = -t175 * t166 + t203;
t96 = t215 * pkin(4) + pkin(9);
t122 = -t173 * t228 + t35 * t67 - t46 * t96;
t121 = t14 * t72 - t39 * t51 + t40 * t52 - t213;
t120 = -qJD(6) * t132 - t3 * t110;
t118 = t120 + t1;
t117 = t71 * t67 - t119;
t101 = t107 * t116;
t97 = -t113 * pkin(4) - pkin(5);
t63 = -qJD(4) * t77 + t178;
t49 = t51 * t175;
t43 = t48 + t172;
t42 = t113 * t65 + t171;
t37 = t217 - t218;
t34 = t153 + t123;
t33 = -t45 + t222;
t28 = t52 * pkin(5) + t51 * pkin(9) + t87;
t25 = t54 * qJD(5) - t113 * t124 + t215 * t63;
t24 = -t53 * qJD(5) + t113 * t63 + t215 * t124;
t22 = t110 * t48 + t112 * t39;
t21 = -t110 * t39 + t112 * t48;
t19 = t110 * t43 + t112 * t42;
t18 = -t110 * t42 + t112 * t43;
t10 = t125 - t209;
t9 = t211 - t226;
t8 = t56 * t157 - t192;
t7 = t156 * t58 - t193;
t6 = -qJD(6) * t27 - t110 * t24 + t112 * t28;
t5 = qJD(6) * t26 + t110 * t28 + t112 * t24;
t4 = (-t29 - t224) * t112 + (-t30 - t159) * t110;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, qJ(2) * t102, 0, 0, 0, 0, 0, 0, 0, t102, t216, t100 * qJD(2) + t85 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t109) * qJD(1), -0.2e1 * t149, 0.2e1 * (t106 - t107) * t176, -t115 * t111, 0.2e1 * t149, -t115 * t114, 0, t111 * t128 + t135 * t182, t114 * t128 - t135 * t183, -0.2e1 * qJD(2) * t154, t155 * qJD(3) + (qJD(1) * t108 + t84) * t219, t138, t45 * t72 + t51 * t67 - t69 * t52 - t206, -t49, t136, -t223, 0, -t25 * t175 + t91 * t46 + t71 * t52 + t87 * t67 + t75 * t72, -t24 * t175 - t91 * t45 - t71 * t51 + t87 * t69 + t75 * t73, -t24 * t67 + t25 * t69 - t53 * t45 - t54 * t46 - t121, t14 * t54 + t40 * t24 - t39 * t25 + t71 * t87 + t75 * t91 + t214, t112 * t142 - t58 * t170 (t110 * t58 + t112 * t56) * t51 + (t193 - t192 + (-t196 + t199) * qJD(6)) * t73, t112 * t137 - t170 * t228 + t141, t140 * t110 + t169 * t56, -t110 * t137 - t169 * t228 + t139, t228 * t52 + t212, t110 * t143 + t228 * t6 + t25 * t56 + t26 * t46 + t53 * t30 + t32 * t73 - t145, t112 * t143 - t228 * t5 + t25 * t58 - t27 * t46 - t53 * t29 - t31 * t73 - t144, t26 * t29 - t27 * t30 - t5 * t56 - t6 * t58 + t132 * t51 + (qJD(6) * t131 - t2 * t110 - t3 * t112) * t73, -t130 * t6 + t17 * t5 + t2 * t27 + t35 * t25 + t3 * t26 + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t116 * qJ(2), 0, 0, 0, 0, 0, 0, 0, -t116, 0, -t100 * qJD(1) - t104, 0, 0, 0, 0, 0, 0, -0.2e1 * t161, 0.2e1 * t162, t101 + t191, -t154 * t84 - t104, 0, 0, 0, 0, 0, 0, -t153 + t123, t45 + t222, t217 + t218, -t39 * t69 - t40 * t67 - t75, 0, 0, 0, 0, 0, 0, t211 + t226, t125 + t209, t232 * t112 + (t30 - t159) * t110, t229 * t110 + t230 * t112 + t35 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, t189 * qJD(1), 0, 0, 0, 0, 0, 0, t186 * t111, t186 * t114, 0 (-t85 + t219) * qJD(1), 0, 0, 0, 0, 0, 0, -qJD(1) * t67 - t49, -qJD(1) * t69 - t223, -t136 - t138, -t71 * qJD(1) + t121, 0, 0, 0, 0, 0, 0, -t72 * t200 + (-t110 * t52 - t112 * t152) * t228 - t140, -t72 * t197 + (t110 * t152 - t112 * t52) * t228 - t142 (t152 * t58 + t139) * t112 + (t152 * t56 + t141) * t110 (t130 * t152 + t144) * t112 + (-t152 * t17 + t145) * t110 - t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t101 - t191, 0, -t168, 0, 0, -t85 * t184 + t94, -t189 * t185, 0, 0, t207, t37, t33, -t207, t34, 0, t41 * qJD(4) - t204 * qJD(5) + (-t175 * t163 - t67 * t184) * pkin(4) + t127, -t69 * t172 + t42 * t175 + (-t175 * t202 - t158) * t113 + t117, t204 * t69 + (-t39 + t42) * t67 + (-t215 * t46 + t113 * t45 + (-t113 * t67 + t215 * t69) * qJD(5)) * pkin(4), t39 * t41 - t40 * t42 + (-t71 * t184 + t215 * t14 - t113 * t15 + (t113 * t40 - t215 * t39) * qJD(5)) * pkin(4), t7, t4, t10, t8, t9, -t208, -t18 * t228 + t97 * t30 + t129 * t56 + (-qJD(6) * t228 * t96 - t15) * t112 + t122 * t110 + t165, -t97 * t29 + (t180 * t96 + t19) * t228 + t129 * t58 + t122 * t112 + t150, t18 * t58 + t19 * t56 + (-t56 * t173 - t30 * t96 + (t58 * t96 + t130) * qJD(6)) * t112 + (t58 * t173 - t29 * t96 - t3 + (t56 * t96 - t17) * qJD(6)) * t110 + t174, t15 * t97 + t130 * t18 - t17 * t19 - t35 * t41 + (-t131 * t113 + t215 * t35) * t202 + t118 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, t37, t33, -t207, t34, 0, t40 * qJD(4) + t127, t39 * t175 + t117 - t221, 0, 0, t7, t4, t10, t8, t9, -t208, t35 * t198 - pkin(5) * t30 - t15 * t112 - t21 * t228 - t40 * t56 + (-t179 * t228 - t200) * pkin(9) + t165, t35 * t195 + pkin(5) * t29 + t22 * t228 - t40 * t58 + (t180 * t228 - t197) * pkin(9) + t150, t21 * t58 + t22 * t56 + (-t193 - t192 + (t196 + t199) * qJD(6)) * pkin(9) + t120 + t174, -t15 * pkin(5) + pkin(9) * t118 + t130 * t21 - t17 * t22 - t35 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, -t56 ^ 2 + t58 ^ 2, t232, -t210, t201 + (-qJD(6) + t228) * t58, t46, -t35 * t58 - t230, t35 * t56 + t229, 0, 0;];
tauc_reg  = t11;
