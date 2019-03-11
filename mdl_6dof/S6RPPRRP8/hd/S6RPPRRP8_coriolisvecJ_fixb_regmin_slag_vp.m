% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:34
% EndTime: 2019-03-09 02:16:37
% DurationCPUTime: 1.87s
% Computational Cost: add. (3540->292), mult. (7752->380), div. (0->0), fcn. (5418->6), ass. (0->137)
t105 = sin(pkin(9));
t106 = cos(pkin(9));
t109 = sin(qJ(4));
t181 = cos(qJ(4));
t119 = -t181 * t105 - t109 * t106;
t189 = t119 * qJD(1);
t193 = qJD(5) - t189;
t108 = sin(qJ(5));
t110 = cos(qJ(5));
t107 = -pkin(1) - qJ(3);
t190 = t107 * qJD(1);
t89 = qJD(2) + t190;
t145 = -pkin(7) * qJD(1) + t89;
t70 = t145 * t105;
t71 = t145 * t106;
t42 = t109 * t71 + t181 * t70;
t38 = qJD(4) * pkin(8) + t42;
t148 = t181 * t106;
t138 = qJD(1) * t148;
t159 = qJD(1) * t105;
t147 = t109 * t159;
t77 = t138 - t147;
t104 = qJD(1) * qJ(2);
t98 = qJD(3) + t104;
t85 = pkin(3) * t159 + t98;
t39 = -pkin(4) * t189 - t77 * pkin(8) + t85;
t12 = t108 * t39 + t110 * t38;
t8 = qJ(6) * t193 + t12;
t195 = t193 * t8;
t58 = t108 * qJD(4) + t110 * t77;
t144 = t58 * t193;
t113 = t189 * qJD(4);
t194 = qJD(5) * qJD(4) + t113;
t118 = -t109 * t105 + t148;
t157 = qJD(5) * t108;
t30 = -t194 * t110 + t77 * t157;
t146 = qJD(4) * t181;
t158 = qJD(4) * t109;
t78 = -t105 * t146 - t106 * t158;
t131 = -t118 * t30 + t78 * t58;
t140 = -qJD(5) * t119 + qJD(1);
t88 = qJD(4) * t147;
t69 = qJD(4) * t138 - t88;
t65 = t110 * t69;
t79 = -t105 * t158 + t106 * t146;
t192 = t119 * t65 + (t108 * t140 - t110 * t79) * t193 - t131;
t191 = -t109 * t70 + t181 * t71;
t174 = -pkin(7) + t107;
t83 = t174 * t105;
t84 = t174 * t106;
t53 = t109 * t83 - t181 * t84;
t160 = t105 ^ 2 + t106 ^ 2;
t188 = t160 * qJD(3);
t156 = qJD(5) * t110;
t114 = t119 * qJD(3);
t21 = qJD(1) * t114 + qJD(4) * t191;
t103 = qJD(1) * qJD(2);
t40 = t69 * pkin(4) - pkin(8) * t113 + t103;
t143 = t108 * t21 - t110 * t40 + t38 * t156 + t39 * t157;
t182 = t69 * pkin(5);
t2 = t143 - t182;
t187 = -t2 + t195;
t37 = -qJD(4) * pkin(4) - t191;
t56 = -t110 * qJD(4) + t108 * t77;
t13 = t56 * pkin(5) - t58 * qJ(6) + t37;
t183 = pkin(8) * t69;
t186 = t13 * t193 - t183;
t184 = t58 ^ 2;
t153 = 0.2e1 * t103;
t180 = t12 * t193;
t179 = t56 * t189;
t178 = t58 * t56;
t177 = t58 * t77;
t176 = t77 * t56;
t175 = t118 * t69;
t127 = pkin(5) * t108 - qJ(6) * t110;
t173 = t108 * qJD(6) - t127 * t193 + t42;
t31 = t108 * t194 + t77 * t156;
t172 = -t108 * t31 - t56 * t156;
t51 = t77 * pkin(4) - pkin(8) * t189;
t171 = t108 * t51 + t110 * t191;
t94 = t105 * pkin(3) + qJ(2);
t50 = -pkin(4) * t119 - pkin(8) * t118 + t94;
t54 = t109 * t84 + t181 * t83;
t170 = t108 * t50 + t110 * t54;
t169 = pkin(8) * qJD(5);
t63 = t108 * t69;
t168 = t110 * t193;
t166 = t30 * t108;
t165 = t69 * qJ(6);
t164 = qJD(5) * t118;
t163 = t79 * qJD(4);
t11 = -t108 * t38 + t110 * t39;
t161 = qJD(6) - t11;
t154 = t193 * t169;
t150 = t118 * t157;
t149 = t118 * t156;
t142 = t108 * t193;
t141 = qJD(1) * t160;
t117 = -qJD(3) * t77 - t70 * t146 - t71 * t158;
t4 = t31 * pkin(5) + t30 * qJ(6) - t58 * qJD(6) - t117;
t139 = -t4 - t154;
t121 = t108 * t40 + t110 * t21 + t39 * t156 - t38 * t157;
t1 = qJD(6) * t193 + t121 + t165;
t137 = -t1 * t119 + t8 * t79;
t7 = -pkin(5) * t193 + t161;
t136 = -t119 * t2 + t7 * t79;
t135 = t118 * t4 + t13 * t78;
t134 = t193 * t7 + t1;
t133 = -t108 * t8 + t110 * t7;
t132 = -t117 * t118 + t37 * t78;
t130 = t119 * t30 + t58 * t79;
t129 = t119 * t31 - t56 * t79;
t128 = -t193 * t78 - t175;
t126 = t156 * t193 - t168 * t189 + t63;
t125 = t65 + (t108 * t189 - t157) * t193;
t124 = t13 * t58 + t143;
t122 = t193 * t37 - t183;
t32 = -t53 * qJD(4) + t114;
t48 = t79 * pkin(4) - t78 * pkin(8) + qJD(2);
t120 = t108 * t48 + t110 * t32 + t50 * t156 - t54 * t157;
t112 = t119 * t63 - t118 * t31 - t78 * t56 + (-t108 * t79 - t140 * t110) * t193;
t33 = t118 * qJD(3) + t54 * qJD(4);
t111 = qJD(1) ^ 2;
t86 = -t110 * pkin(5) - t108 * qJ(6) - pkin(4);
t72 = t78 * qJD(4);
t24 = pkin(5) * t58 + qJ(6) * t56;
t20 = t118 * t127 + t53;
t16 = pkin(5) * t119 + t108 * t54 - t110 * t50;
t15 = -qJ(6) * t119 + t170;
t14 = t193 * t56 - t30;
t10 = -t77 * pkin(5) + t108 * t191 - t110 * t51;
t9 = t77 * qJ(6) + t171;
t6 = (pkin(5) * t78 + qJ(6) * t164) * t108 + (-qJ(6) * t78 - (-pkin(5) * qJD(5) + qJD(6)) * t118) * t110 + t33;
t5 = -t79 * pkin(5) + t170 * qJD(5) + t108 * t32 - t110 * t48;
t3 = t79 * qJ(6) - qJD(6) * t119 + t120;
t17 = [0, 0, 0, 0, t153, qJ(2) * t153, t105 * t153, t106 * t153, 0.2e1 * qJD(3) * t141 (t98 + t104) * qJD(2) + (-t89 - t190) * t188, t113 * t118 + t77 * t78, t113 * t119 + t189 * t78 - t77 * t79 - t175, t72, -t163, 0, -0.2e1 * qJD(2) * t189 - t33 * qJD(4) + t94 * t69 + t85 * t79, qJD(2) * t77 - t32 * qJD(4) + t103 * t118 + t113 * t94 + t85 * t78, t110 * t131 - t150 * t58 (-t108 * t58 - t110 * t56) * t78 - (-t166 + t110 * t31 + (-t108 * t56 + t110 * t58) * qJD(5)) * t118, -t110 * t128 - t150 * t193 + t130, t108 * t128 - t149 * t193 + t129, -t119 * t69 + t193 * t79, t143 * t119 + t11 * t79 + t33 * t56 + t53 * t31 + ((-qJD(5) * t54 + t48) * t193 + t50 * t69 + t37 * t164) * t110 + ((-qJD(5) * t50 - t32) * t193 - t54 * t69 + t132) * t108, t110 * t132 + t119 * t121 - t12 * t79 - t120 * t193 - t150 * t37 - t170 * t69 - t53 * t30 + t33 * t58, t108 * t135 + t13 * t149 - t16 * t69 - t193 * t5 + t20 * t31 + t6 * t56 - t136, -t15 * t31 - t16 * t30 - t3 * t56 + t5 * t58 + t133 * t78 - (t1 * t108 - t2 * t110 + (t108 * t7 + t110 * t8) * qJD(5)) * t118, -t110 * t135 + t13 * t150 + t15 * t69 + t193 * t3 + t20 * t30 - t6 * t58 + t137, t1 * t15 + t13 * t6 + t16 * t2 + t20 * t4 + t3 * t8 + t5 * t7; 0, 0, 0, 0, -t111, -t111 * qJ(2), -t111 * t105, -t111 * t106, 0 (-t98 - t188) * qJD(1), 0, 0, 0, 0, 0, qJD(1) * t189 + t72, -qJD(1) * t77 - t163, 0, 0, 0, 0, 0, t112, t192, t112 (t140 * t58 + t129) * t110 + (t140 * t56 + t130) * t108, -t192 (t140 * t7 + t137) * t110 + (-t140 * t8 + t136) * t108 - t135; 0, 0, 0, 0, 0, 0, 0, 0, -t160 * t111, t89 * t141 + t103, 0, 0, 0, 0, 0, -t88 + (t138 + t77) * qJD(4), 0.2e1 * t113, 0, 0, 0, 0, 0, t125 - t176, -t168 * t193 - t177 - t63, -t142 * t193 - t176 + t65 (t30 + t179) * t110 + t108 * t144 + t172, t126 + t177, t134 * t108 + t187 * t110 - t13 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77 * t189, -t189 ^ 2 + t77 ^ 2, 0, t88 + (-t138 + t77) * qJD(4), 0, t42 * qJD(4) - t85 * t77 + t117 (-qJD(3) - t85) * t189, t110 * t144 - t166 (-t30 + t179) * t110 - t58 * t142 + t172, t126 - t177, t125 + t176, -t193 * t77, -pkin(4) * t31 - t11 * t77 - t42 * t56 + (t117 + (-t51 - t169) * t193) * t110 + (t191 * t193 + t122) * t108, pkin(4) * t30 + t171 * t193 + t12 * t77 - t42 * t58 + (-t117 + t154) * t108 + t122 * t110, t10 * t193 + t186 * t108 + t139 * t110 - t173 * t56 + t86 * t31 + t7 * t77, -t10 * t58 + t9 * t56 + ((qJD(5) * t58 - t31) * pkin(8) + t134) * t110 + ((qJD(5) * t56 - t30) * pkin(8) - t187) * t108, t139 * t108 - t186 * t110 + t173 * t58 - t193 * t9 + t86 * t30 - t8 * t77, -t7 * t10 + t4 * t86 - t8 * t9 - t173 * t13 + (qJD(5) * t133 + t1 * t110 + t2 * t108) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, -t56 ^ 2 + t184, t14, -t31 + t144, t69, -t37 * t58 - t143 + t180, t11 * t193 + t37 * t56 - t121, -t24 * t56 - t124 + t180 + 0.2e1 * t182, pkin(5) * t30 - qJ(6) * t31 + (-t12 + t8) * t58 + (t7 - t161) * t56, 0.2e1 * t165 - t13 * t56 + t24 * t58 + (0.2e1 * qJD(6) - t11) * t193 + t121, -pkin(5) * t2 + qJ(6) * t1 - t12 * t7 - t13 * t24 + t161 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 + t178, t14, -t193 ^ 2 - t184, t124 - t182 - t195;];
tauc_reg  = t17;
