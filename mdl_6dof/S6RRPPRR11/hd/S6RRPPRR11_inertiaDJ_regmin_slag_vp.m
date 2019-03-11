% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:42:44
% EndTime: 2019-03-09 09:42:52
% DurationCPUTime: 2.37s
% Computational Cost: add. (2982->273), mult. (7820->518), div. (0->0), fcn. (7717->10), ass. (0->149)
t113 = sin(pkin(6));
t119 = sin(qJ(2));
t116 = -pkin(2) - qJ(4);
t121 = cos(qJ(2));
t179 = qJ(3) * t119;
t135 = -t116 * t121 + t179;
t197 = t113 * (qJD(2) * t135 + qJD(4) * t119);
t115 = cos(pkin(6));
t186 = pkin(1) * t115;
t101 = t119 * t186;
t176 = t113 * t121;
t196 = pkin(8) * t176 + t101;
t112 = sin(pkin(11));
t114 = cos(pkin(11));
t118 = sin(qJ(5));
t187 = cos(qJ(5));
t195 = -t118 * t112 + t187 * t114;
t74 = t115 * t112 + t114 * t176;
t75 = -t112 * t176 + t115 * t114;
t129 = -t118 * t75 - t187 * t74;
t169 = qJD(2) * t119;
t153 = t113 * t169;
t81 = t187 * t112 + t118 * t114;
t32 = qJD(5) * t129 + t81 * t153;
t147 = t114 * t153;
t148 = t112 * t153;
t49 = -t118 * t74 + t187 * t75;
t33 = t49 * qJD(5) + t118 * t148 - t187 * t147;
t194 = t33 * pkin(5) - t32 * pkin(10);
t120 = cos(qJ(6));
t111 = t120 ^ 2;
t117 = sin(qJ(6));
t172 = t117 ^ 2 - t111;
t149 = qJD(6) * t172;
t87 = (t112 ^ 2 + t114 ^ 2) * qJD(4);
t177 = t113 * t119;
t156 = -pkin(1) * t121 - pkin(2);
t98 = pkin(8) * t177;
t57 = pkin(3) * t177 + t98 + (-qJ(4) + t156) * t115;
t63 = (-pkin(1) - t135) * t113;
t34 = -t112 * t63 + t114 * t57;
t25 = pkin(4) * t177 - t75 * pkin(9) + t34;
t35 = t112 * t57 + t114 * t63;
t27 = -t74 * pkin(9) + t35;
t130 = t118 * t25 + t187 * t27;
t170 = qJD(2) * t113;
t178 = t112 * t119;
t167 = qJD(3) * t119;
t91 = pkin(2) * t153;
t45 = t91 + (-t167 - qJD(4) * t121 + (-qJ(3) * t121 + qJ(4) * t119) * qJD(2)) * t113;
t190 = pkin(3) + pkin(8);
t59 = -t115 * qJD(4) + (t190 * t176 + t101) * qJD(2);
t28 = -t112 * t45 + t114 * t59;
t20 = (pkin(4) * t121 - pkin(9) * t178) * t170 + t28;
t29 = t112 * t59 + t114 * t45;
t24 = pkin(9) * t147 + t29;
t6 = -qJD(5) * t130 - t118 * t24 + t187 * t20;
t79 = t195 ^ 2;
t192 = 0.2e1 * t113;
t191 = 0.2e1 * qJD(3);
t151 = qJD(5) * t187;
t165 = qJD(5) * t118;
t77 = -t112 * t151 - t114 * t165;
t185 = t195 * t77;
t78 = -t112 * t165 + t114 * t151;
t184 = t81 * t78;
t183 = -pkin(9) + t116;
t124 = -t117 * t49 + t120 * t177;
t168 = qJD(2) * t121;
t94 = t113 * t168;
t14 = qJD(6) * t124 + t117 * t94 + t120 * t32;
t182 = t14 * t117;
t73 = t196 * qJD(2);
t181 = t73 * t115;
t104 = t115 * qJD(3);
t92 = t168 * t186;
t180 = t104 + t92;
t175 = t117 * t120;
t102 = t112 * pkin(4) + qJ(3);
t171 = qJ(3) * qJD(3);
t164 = qJD(6) * t117;
t163 = qJD(6) * t120;
t162 = -0.2e1 * pkin(5) * qJD(6);
t161 = t81 ^ 2 + t79;
t160 = t81 * t164;
t159 = t81 * t163;
t158 = t195 * t164;
t157 = t195 * t163;
t67 = -t115 * qJ(3) - t196;
t108 = t113 ^ 2;
t154 = t108 * t168;
t152 = t117 * t163;
t150 = -0.4e1 * t195 * t175;
t64 = pkin(3) * t176 - t67;
t146 = -pkin(5) * t77 - pkin(10) * t78;
t145 = -pkin(5) * t195 - pkin(10) * t81;
t4 = -pkin(5) * t94 - t6;
t132 = -t118 * t27 + t187 * t25;
t9 = -pkin(5) * t177 - t132;
t144 = t195 * t4 + t77 * t9;
t37 = t117 * t177 + t120 * t49;
t143 = t14 * t195 + t77 * t37;
t142 = t129 * t77 - t195 * t33;
t141 = -t129 * t78 + t33 * t81;
t84 = t183 * t112;
t85 = t183 * t114;
t62 = t118 * t85 + t187 * t84;
t40 = qJD(4) * t195 + t62 * qJD(5);
t61 = t118 * t84 - t187 * t85;
t140 = t195 * t40 + t61 * t77;
t139 = t195 * t78 + t77 * t81;
t138 = -pkin(2) * t121 - t179;
t10 = pkin(10) * t177 + t130;
t42 = t74 * pkin(4) + t64;
t17 = -pkin(5) * t129 - t49 * pkin(10) + t42;
t8 = t120 * t10 + t117 * t17;
t137 = t29 * t112 + t28 * t114;
t136 = -t117 * t37 + t120 * t124;
t52 = t81 * pkin(5) - pkin(10) * t195 + t102;
t31 = t117 * t52 + t120 * t62;
t134 = (-pkin(4) * t114 - t190) * t119;
t72 = pkin(8) * t153 - t92;
t133 = t78 * pkin(5) - t77 * pkin(10) + qJD(3);
t128 = t117 * t33 - t129 * t163;
t127 = -t120 * t33 - t129 * t164;
t126 = -t117 * t77 - t157;
t125 = t120 * t77 - t158;
t51 = t117 * t78 + t159;
t5 = -t118 * t20 - t25 * t151 + t27 * t165 - t187 * t24;
t123 = -0.2e1 * t184 - 0.2e1 * t185;
t122 = (-t119 * t78 - t81 * t168) * t113;
t44 = t134 * t170 + t180;
t86 = 0.2e1 * t119 * t154;
t70 = t156 * t115 + t98;
t68 = (-pkin(1) + t138) * t113;
t66 = -t104 + t72;
t65 = t91 + (-qJ(3) * t168 - t167) * t113;
t58 = -t190 * t153 + t180;
t50 = -t120 * t78 + t160;
t43 = (t119 * t77 + t168 * t195) * t113;
t39 = t81 * qJD(4) - t85 * t151 + t84 * t165;
t30 = -t117 * t62 + t120 * t52;
t15 = qJD(6) * t37 + t117 * t32 - t120 * t94;
t13 = -qJD(6) * t31 + t117 * t39 + t120 * t133;
t12 = -t117 * t133 + t120 * t39 - t52 * t163 + t62 * t164;
t7 = -t117 * t10 + t120 * t17;
t2 = t117 * t5 + t120 * (t180 + t194) - t8 * qJD(6) + (-t117 * pkin(10) * t121 + t120 * t134) * t170;
t1 = t10 * t164 - t117 * (t44 + t194) - t17 * t163 - t120 * (pkin(10) * t94 - t5);
t3 = [0, 0, 0, t86, 0.2e1 * (-t119 ^ 2 + t121 ^ 2) * t108 * qJD(2), 0.2e1 * t115 * t94, -0.2e1 * t115 * t153, 0, -0.2e1 * t108 * pkin(1) * t169 - 0.2e1 * t181, -0.2e1 * pkin(1) * t154 + 0.2e1 * t72 * t115 (t119 * t73 - t121 * t66 + (t119 * t67 + t121 * t70) * qJD(2)) * t192, 0.2e1 * t181 + 0.2e1 * (t121 * t65 - t68 * t169) * t113, -0.2e1 * t66 * t115 + 0.2e1 * (-t119 * t65 - t68 * t168) * t113, 0.2e1 * t68 * t65 + 0.2e1 * t67 * t66 + 0.2e1 * t70 * t73, 0.2e1 * t58 * t74 + 0.2e1 * (t119 * t28 + (-t114 * t119 * t64 + t121 * t34) * qJD(2)) * t113, 0.2e1 * t58 * t75 + 0.2e1 * (-t119 * t29 + (-t121 * t35 + t64 * t178) * qJD(2)) * t113, -0.2e1 * t28 * t75 - 0.2e1 * t29 * t74 + 0.2e1 * (-t112 * t34 + t114 * t35) * t153, 0.2e1 * t34 * t28 + 0.2e1 * t35 * t29 + 0.2e1 * t64 * t58, 0.2e1 * t49 * t32, 0.2e1 * t129 * t32 - 0.2e1 * t49 * t33 (t119 * t32 + t49 * t168) * t192 (-t119 * t33 + t129 * t168) * t192, t86, 0.2e1 * t42 * t33 - 0.2e1 * t44 * t129 + 0.2e1 * (t6 * t119 + t132 * t168) * t113, 0.2e1 * t42 * t32 + 0.2e1 * t44 * t49 + 0.2e1 * (t5 * t119 - t130 * t168) * t113, 0.2e1 * t37 * t14, 0.2e1 * t124 * t14 - 0.2e1 * t37 * t15, -0.2e1 * t129 * t14 + 0.2e1 * t37 * t33, 0.2e1 * t124 * t33 + 0.2e1 * t129 * t15, -0.2e1 * t129 * t33, -0.2e1 * t124 * t4 - 0.2e1 * t129 * t2 + 0.2e1 * t9 * t15 + 0.2e1 * t7 * t33, -0.2e1 * t1 * t129 + 0.2e1 * t9 * t14 - 0.2e1 * t8 * t33 + 0.2e1 * t4 * t37; 0, 0, 0, 0, 0, t94, -t153, 0, -t73, t72 (t138 * qJD(2) + qJD(3) * t121) * t113, t73, 0.2e1 * t104 - t72, -t73 * pkin(2) - t66 * qJ(3) - t67 * qJD(3), qJD(3) * t74 + t58 * t112 - t114 * t197, qJD(3) * t75 + t112 * t197 + t58 * t114 (t112 * t74 + t114 * t75) * qJD(4) - t137, t58 * qJ(3) + t64 * qJD(3) + t137 * t116 + (-t112 * t35 - t114 * t34) * qJD(4), t195 * t32 + t49 * t77, -t32 * t81 - t49 * t78 + t142, t43, t122, 0, -qJD(3) * t129 + t102 * t33 + t42 * t78 + t44 * t81 + (-t119 * t40 - t61 * t168) * t113, qJD(3) * t49 + t102 * t32 + t42 * t77 + t44 * t195 + (t119 * t39 - t62 * t168) * t113, t143 * t120 - t37 * t158, t136 * t77 - (t182 + t120 * t15 + (t117 * t124 + t120 * t37) * qJD(6)) * t195, -t142 * t120 + t129 * t158 + t14 * t81 + t37 * t78, t142 * t117 + t124 * t78 + t129 * t157 - t15 * t81, t141, t144 * t117 - t124 * t40 - t129 * t13 + t61 * t15 + t9 * t157 + t2 * t81 + t30 * t33 + t7 * t78, t1 * t81 - t12 * t129 + t144 * t120 + t61 * t14 - t9 * t158 - t31 * t33 + t40 * t37 - t8 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, 0.2e1 * t171, t112 * t191, t114 * t191, 0.2e1 * t87, -0.2e1 * t116 * t87 + 0.2e1 * t171, 0.2e1 * t185, -0.2e1 * t139, 0, 0, 0, 0.2e1 * qJD(3) * t81 + 0.2e1 * t102 * t78, 0.2e1 * qJD(3) * t195 + 0.2e1 * t102 * t77, 0.2e1 * t111 * t185 - 0.2e1 * t79 * t152, 0.2e1 * t79 * t149 + t77 * t150, 0.2e1 * t139 * t120 - 0.2e1 * t81 * t158, -0.2e1 * t139 * t117 - 0.2e1 * t81 * t157, 0.2e1 * t184, 0.2e1 * t140 * t117 + 0.2e1 * t13 * t81 + 0.2e1 * t61 * t157 + 0.2e1 * t30 * t78, 0.2e1 * t12 * t81 + 0.2e1 * t140 * t120 - 0.2e1 * t61 * t158 - 0.2e1 * t31 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, 0, t73, t114 * t94, -t112 * t94, 0, t137, 0, 0, 0, 0, 0, t43, t122, 0, 0, 0, 0, 0, -t141 * t117 + t124 * t77 + t129 * t159 - t15 * t195, -t141 * t120 - t129 * t160 - t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117 * t123 - t161 * t163, t120 * t123 + t161 * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t148, 0, t58, 0, 0, 0, 0, 0, t33, t32, 0, 0, 0, 0, 0, -t127, -t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, t78, t77, 0, 0, 0, 0, 0, -t50, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t33, t94, t6, t5, t37 * t163 + t182, t136 * qJD(6) - t117 * t15 + t14 * t120, t128, -t127, 0, -pkin(5) * t15 - pkin(10) * t128 - t4 * t120 + t9 * t164, -pkin(5) * t14 + pkin(10) * t127 + t4 * t117 + t9 * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t78, 0, -t40, t39, -t149 * t195 + t77 * t175, qJD(6) * t150 - t172 * t77, t51, -t50, 0, -t40 * t120 + t146 * t117 + (t117 * t61 + t145 * t120) * qJD(6), t40 * t117 + t146 * t120 + (-t145 * t117 + t120 * t61) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t78, 0, 0, 0, 0, 0, t125, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t152, -0.2e1 * t149, 0, 0, 0, t117 * t162, t120 * t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, t33, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t126, t78, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, -t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, -t164, 0, -pkin(10) * t163, pkin(10) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
