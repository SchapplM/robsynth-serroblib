% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:06:31
% EndTime: 2019-03-08 21:06:35
% DurationCPUTime: 1.78s
% Computational Cost: add. (2167->271), mult. (5768->385), div. (0->0), fcn. (4383->10), ass. (0->160)
t108 = sin(pkin(11));
t110 = cos(pkin(11));
t113 = sin(qJ(3));
t116 = cos(qJ(3));
t87 = t108 * t116 + t110 * t113;
t82 = t87 * qJD(2);
t207 = qJD(6) + t82;
t112 = sin(qJ(6));
t115 = cos(qJ(6));
t166 = qJD(2) * t113;
t175 = t110 * t116;
t79 = -qJD(2) * t175 + t108 * t166;
t60 = qJD(3) * t112 - t115 * t79;
t208 = t207 * t60;
t147 = t112 * t207;
t160 = qJD(2) * qJD(3);
t150 = t116 * t160;
t151 = t113 * t160;
t71 = -t108 * t151 + t110 * t150;
t133 = t115 * t71 - t147 * t207;
t114 = sin(qJ(2));
t109 = sin(pkin(6));
t168 = qJD(1) * t109;
t156 = t114 * t168;
t164 = qJD(3) * t113;
t206 = pkin(3) * t164 - t156;
t205 = -0.2e1 * qJD(3);
t117 = cos(qJ(2));
t153 = t117 * t168;
t187 = -qJ(4) - pkin(8);
t149 = qJD(3) * t187;
t75 = t116 * qJD(4) + t113 * t149;
t76 = -t113 * qJD(4) + t116 * t149;
t184 = t108 * t75 - t110 * t76 - t87 * t153;
t200 = -t108 * t113 + t175;
t183 = t108 * t76 + t110 * t75 - t200 * t153;
t78 = t82 ^ 2;
t203 = -t79 ^ 2 - t78;
t90 = qJD(2) * pkin(8) + t156;
t146 = qJ(4) * qJD(2) + t90;
t111 = cos(pkin(6));
t167 = qJD(1) * t111;
t152 = t113 * t167;
t55 = t146 * t116 + t152;
t47 = t108 * t55;
t98 = t116 * t167;
t54 = -t146 * t113 + t98;
t24 = t110 * t54 - t47;
t172 = -qJD(5) + t24;
t202 = -qJD(6) + t207;
t163 = qJD(3) * t116;
t84 = -t108 * t164 + t110 * t163;
t201 = qJ(5) * t84 + qJD(5) * t87 - t206;
t181 = t110 * t55;
t51 = qJD(3) * pkin(3) + t54;
t20 = t108 * t51 + t181;
t17 = -qJD(3) * qJ(5) - t20;
t195 = pkin(5) * t79;
t11 = -t17 - t195;
t22 = t108 * t54 + t181;
t103 = -pkin(3) * t110 - pkin(4);
t99 = -pkin(9) + t103;
t199 = t99 * t71 + (t11 - t22 + t195) * t207;
t197 = pkin(4) + pkin(9);
t81 = t87 * qJD(3);
t70 = qJD(2) * t81;
t196 = pkin(4) * t70;
t194 = pkin(5) * t82;
t178 = t109 * t114;
t128 = t111 * t116 - t113 * t178;
t85 = t111 * t113 + t116 * t178;
t43 = t108 * t85 - t110 * t128;
t136 = qJD(4) + t153;
t33 = (-t113 * t90 + t98) * qJD(3) + (-qJ(4) * t164 + t136 * t116) * qJD(2);
t34 = (-t116 * t90 - t152) * qJD(3) + (-qJ(4) * t163 - t136 * t113) * qJD(2);
t7 = t108 * t33 - t110 * t34;
t193 = t43 * t7;
t92 = t187 * t113;
t93 = t187 * t116;
t58 = -t108 * t93 - t110 * t92;
t192 = t7 * t58;
t191 = t11 * t200;
t157 = -pkin(3) * t116 - pkin(2);
t135 = -qJ(5) * t87 + t157;
t35 = -t197 * t200 + t135;
t190 = t35 * t71;
t62 = qJD(3) * t115 + t112 * t79;
t189 = t62 * t79;
t188 = t79 * t60;
t186 = pkin(5) * t84 + t184;
t185 = -pkin(5) * t81 + t183;
t8 = t108 * t34 + t110 * t33;
t182 = qJD(2) * pkin(2);
t161 = qJD(6) * t115;
t162 = qJD(6) * t112;
t36 = -qJD(3) * t162 + t112 * t70 + t79 * t161;
t180 = t36 * t115;
t165 = qJD(2) * t114;
t155 = t109 * t165;
t77 = pkin(3) * t151 + qJD(1) * t155;
t177 = t109 * t117;
t119 = qJD(2) ^ 2;
t176 = t109 * t119;
t118 = qJD(3) ^ 2;
t174 = t118 * t113;
t173 = t118 * t116;
t171 = t194 - t172;
t169 = t113 ^ 2 - t116 ^ 2;
t159 = t200 * t161;
t158 = t114 * t176;
t154 = qJD(2) * t177;
t19 = t110 * t51 - t47;
t148 = pkin(3) * t166 + t79 * qJ(5);
t145 = t113 * t154;
t144 = t116 * t154;
t143 = -qJ(5) * t71 + t77;
t142 = qJD(5) - t19;
t141 = -t197 * t81 + t201;
t140 = -pkin(4) * t81 + t201;
t137 = -t200 * t71 + t207 * t81;
t6 = -qJD(3) * qJD(5) - t8;
t10 = -t197 * qJD(3) + t142 + t194;
t72 = t157 * qJD(2) + qJD(4) - t153;
t121 = -qJ(5) * t82 + t72;
t18 = t197 * t79 + t121;
t1 = t10 * t115 - t112 * t18;
t2 = t10 * t112 + t115 * t18;
t59 = t108 * t92 - t110 * t93;
t32 = pkin(4) * t79 + t121;
t134 = t32 * t82 + t7;
t132 = -t112 * t43 + t115 * t177;
t131 = t112 * t177 + t115 * t43;
t3 = -pkin(5) * t70 - t6;
t129 = t3 + (-qJD(6) * t99 + t197 * t82 + t148) * t207;
t127 = t182 * qJD(2);
t41 = pkin(5) * t87 + t58;
t126 = -t11 * t81 + t200 * t3 + t41 * t71;
t125 = -qJD(5) * t82 + t143;
t124 = -t115 * t207 ^ 2 - t112 * t71;
t123 = t182 * t205;
t52 = t128 * qJD(3) + t144;
t53 = -t85 * qJD(3) - t145;
t21 = t108 * t52 - t110 * t53;
t23 = t108 * t53 + t110 * t52;
t44 = t108 * t128 + t110 * t85;
t122 = t21 * t82 - t23 * t79 + t43 * t71 - t44 * t70;
t120 = -t183 * t79 + t184 * t82 + t58 * t71 - t59 * t70 + t7 * t87;
t101 = pkin(3) * t108 + qJ(5);
t73 = qJD(3) * t79;
t64 = t115 * t70;
t45 = -pkin(4) * t200 + t135;
t42 = pkin(5) * t200 + t59;
t38 = pkin(4) * t82 + t148;
t37 = qJD(6) * t62 - t64;
t16 = -qJD(3) * pkin(4) + t142;
t14 = t125 + t196;
t9 = t197 * t70 + t125;
t5 = pkin(5) * t71 + t7;
t4 = t115 * t5;
t12 = [0, 0, -t158, -t117 * t176, 0, 0, 0, 0, 0, -t116 * t158 + (t53 - t145) * qJD(3), t113 * t158 + (-t52 - t144) * qJD(3), t122, -t19 * t21 + t20 * t23 + t193 + t44 * t8 + (-t117 * t77 + t72 * t165) * t109, t122, t21 * qJD(3) + (t117 * t70 - t79 * t165) * t109, t23 * qJD(3) + (t117 * t71 - t82 * t165) * t109, t16 * t21 - t17 * t23 + t193 - t44 * t6 + (-t117 * t14 + t32 * t165) * t109, 0, 0, 0, 0, 0 (qJD(6) * t132 - t112 * t155 + t115 * t21) * t207 + t131 * t71 + t23 * t60 + t44 * t37 -(qJD(6) * t131 + t112 * t21 + t115 * t155) * t207 + t132 * t71 + t23 * t62 + t44 * t36; 0, 0, 0, 0, 0.2e1 * t113 * t150, -0.2e1 * t169 * t160, t173, -t174, 0, -pkin(8) * t173 + t113 * t123, pkin(8) * t174 + t116 * t123, -t19 * t84 - t20 * t81 + t200 * t8 + t120, t77 * t157 + t183 * t20 - t184 * t19 + t206 * t72 + t8 * t59 + t192, t16 * t84 + t17 * t81 - t200 * t6 + t120, t184 * qJD(3) + t14 * t200 + t140 * t79 - t32 * t81 - t45 * t70, t183 * qJD(3) - t14 * t87 + t140 * t82 - t32 * t84 - t45 * t71, t14 * t45 - t140 * t32 + t184 * t16 - t183 * t17 - t59 * t6 + t192, -t62 * t159 + (-t200 * t36 + t62 * t81) * t112 (-t112 * t60 + t115 * t62) * t81 - (-t112 * t37 + t180 + (-t112 * t62 - t115 * t60) * qJD(6)) * t200, t112 * t137 - t159 * t207 + t36 * t87 + t62 * t84, t162 * t200 * t207 + t115 * t137 - t37 * t87 - t60 * t84, t207 * t84 + t71 * t87, t1 * t84 + t42 * t37 + t4 * t87 + t185 * t60 + (t141 * t207 - t9 * t87 - t190) * t112 + (t186 * t207 + t126) * t115 + ((-t112 * t41 - t115 * t35) * t207 - t2 * t87 - t112 * t191) * qJD(6), -t2 * t84 + t42 * t36 + t185 * t62 + (-t190 - (qJD(6) * t10 + t9) * t87 - qJD(6) * t191 + (-qJD(6) * t41 + t141) * t207) * t115 + (-(-qJD(6) * t18 + t5) * t87 + (qJD(6) * t35 - t186) * t207 - t126) * t112; 0, 0, 0, 0, -t113 * t119 * t116, t169 * t119, 0, 0, 0, t113 * t127, t116 * t127 (t20 - t22) * t82 + (-t19 + t24) * t79 + (-t108 * t70 - t110 * t71) * pkin(3), t19 * t22 - t20 * t24 + (t108 * t8 - t110 * t7 - t72 * t166) * pkin(3), -t101 * t70 + t103 * t71 + (-t17 - t22) * t82 + (t16 + t172) * t79, -qJD(3) * t22 + t38 * t79 + t134, -t32 * t79 + t38 * t82 + (0.2e1 * qJD(5) - t24) * qJD(3) + t8, -t101 * t6 + t103 * t7 - t16 * t22 + t172 * t17 - t32 * t38, -t147 * t62 + t180 (-t207 * t62 - t37) * t115 + (-t36 + t208) * t112, t133 + t189, t124 - t188, t207 * t79, t1 * t79 + t101 * t37 + t129 * t112 + t199 * t115 + t171 * t60, t101 * t36 - t199 * t112 + t129 * t115 + t171 * t62 - t2 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, t19 * t82 + t20 * t79 + t77, t203, t82 * t205, -t71 + t73, t196 - t17 * t79 + (-qJD(5) - t16) * t82 + t143, 0, 0, 0, 0, 0, t124 + t188, t189 - t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 + t73, -t82 * t79, -t78 - t118, qJD(3) * t17 + t134, 0, 0, 0, 0, 0, -qJD(3) * t60 + t133, -qJD(3) * t62 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60, -t60 ^ 2 + t62 ^ 2, t36 + t208, t202 * t62 + t64, t71, -t11 * t62 - t112 * t9 + t202 * t2 + t4, t202 * t1 + t11 * t60 - t112 * t5 - t115 * t9;];
tauc_reg  = t12;
