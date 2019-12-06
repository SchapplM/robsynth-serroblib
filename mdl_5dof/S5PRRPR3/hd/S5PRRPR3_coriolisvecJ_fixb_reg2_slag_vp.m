% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:41
% EndTime: 2019-12-05 16:19:46
% DurationCPUTime: 1.35s
% Computational Cost: add. (2199->208), mult. (5605->292), div. (0->0), fcn. (4030->6), ass. (0->136)
t122 = sin(qJ(5));
t120 = sin(pkin(9));
t121 = cos(pkin(9));
t123 = sin(qJ(3));
t124 = cos(qJ(3));
t100 = t120 * t124 + t121 * t123;
t131 = qJD(2) * t100;
t172 = cos(qJ(5));
t151 = qJD(2) * t124;
t143 = t121 * t151;
t152 = qJD(2) * t123;
t92 = t120 * t152 - t143;
t134 = t122 * t92 - t131 * t172;
t80 = qJD(3) * t131;
t147 = qJD(2) * qJD(3);
t141 = t123 * t147;
t109 = t120 * t141;
t140 = t124 * t147;
t81 = t121 * t140 - t109;
t129 = qJD(5) * t134 - t122 * t81 - t172 * t80;
t117 = qJD(3) + qJD(5);
t160 = t134 * t117;
t185 = t129 - t160;
t142 = qJD(5) * t172;
t150 = qJD(5) * t122;
t133 = -t122 * t80 - t131 * t150 - t142 * t92 + t172 * t81;
t41 = -t122 * t131 - t172 * t92;
t159 = t41 * t117;
t184 = t133 - t159;
t168 = t134 ^ 2;
t169 = t41 ^ 2;
t183 = t168 - t169;
t167 = t41 * t134;
t148 = qJD(1) * qJD(3);
t114 = t124 * t148;
t165 = -qJ(4) - pkin(6);
t138 = qJD(3) * t165;
t87 = t124 * qJD(4) + t123 * t138;
t52 = qJD(2) * t87 + t114;
t89 = -t123 * qJD(4) + t124 * t138;
t53 = qJD(2) * t89 - t123 * t148;
t24 = -t120 * t52 + t121 * t53;
t16 = -pkin(7) * t81 + t24;
t25 = t120 * t53 + t121 * t52;
t17 = -pkin(7) * t80 + t25;
t176 = t131 * pkin(7);
t108 = t165 * t124;
t149 = t123 * qJD(1);
t88 = -qJD(2) * t108 + t149;
t65 = t120 * t88;
t163 = qJD(3) * pkin(3);
t107 = t165 * t123;
t116 = t124 * qJD(1);
t86 = qJD(2) * t107 + t116;
t74 = t86 + t163;
t30 = t121 * t74 - t65;
t22 = qJD(3) * pkin(4) - t176 + t30;
t177 = t92 * pkin(7);
t162 = t121 * t88;
t31 = t120 * t74 + t162;
t23 = t31 - t177;
t6 = t122 * t22 + t172 * t23;
t2 = -qJD(5) * t6 - t122 * t17 + t16 * t172;
t115 = -pkin(3) * t124 - pkin(2);
t153 = qJD(2) * t115;
t106 = qJD(4) + t153;
t50 = t92 * pkin(4) + t106;
t182 = t134 * t50 + t2;
t128 = -t122 * t16 - t142 * t22 + t150 * t23 - t17 * t172;
t181 = -t50 * t41 + t128;
t180 = -0.2e1 * t147;
t179 = 0.2e1 * t131;
t178 = t131 ^ 2;
t125 = qJD(3) ^ 2;
t132 = t100 * qJD(3);
t99 = t120 * t123 - t121 * t124;
t96 = t99 * qJD(3);
t19 = t100 * t150 + t122 * t132 + t142 * t99 + t172 * t96;
t44 = t100 * t172 - t122 * t99;
t175 = t129 * t44 - t19 * t41;
t32 = -t120 * t86 - t162;
t26 = t32 + t177;
t34 = t121 * t86 - t65;
t28 = t34 - t176;
t113 = pkin(3) * t121 + pkin(4);
t171 = pkin(3) * t120;
t85 = t113 * t122 + t171 * t172;
t174 = qJD(5) * t85 - t122 * t28 + t172 * t26;
t84 = t113 * t172 - t122 * t171;
t173 = qJD(5) * t84 - t122 * t26 - t172 * t28;
t170 = t123 * pkin(3);
t166 = t131 * t92;
t164 = -t100 * t80 + t92 * t96;
t35 = t120 * t89 + t121 * t87;
t161 = t19 * t117;
t126 = qJD(2) ^ 2;
t158 = t124 * t126;
t157 = t125 * t123;
t156 = t125 * t124;
t155 = t96 * qJD(3);
t55 = t107 * t120 - t108 * t121;
t154 = t123 ^ 2 - t124 ^ 2;
t146 = pkin(3) * t152;
t145 = pkin(6) * t152;
t144 = t123 * t158;
t111 = pkin(3) * t141;
t51 = pkin(4) * t80 + t111;
t33 = -t120 * t87 + t121 * t89;
t137 = pkin(2) * t180;
t54 = t107 * t121 + t108 * t120;
t136 = t123 * t140;
t20 = qJD(5) * t44 - t122 * t96 + t132 * t172;
t43 = t100 * t122 + t172 * t99;
t135 = t133 * t43 - t134 * t20;
t36 = -pkin(7) * t100 + t54;
t37 = -pkin(7) * t99 + t55;
t11 = -t122 * t37 + t172 * t36;
t12 = t122 * t36 + t172 * t37;
t105 = pkin(6) * t151 + t149;
t130 = -t131 * t132 - t81 * t99;
t104 = t116 - t145;
t97 = -pkin(6) * t141 + t114;
t98 = t105 * qJD(3);
t127 = t98 * t123 + t97 * t124 + (-t104 * t124 - t105 * t123) * qJD(3);
t90 = t92 ^ 2;
t82 = t100 * t125;
t64 = pkin(4) * t99 + t115;
t57 = (pkin(4) * t100 + t170) * qJD(3);
t56 = pkin(4) * t131 + t146;
t29 = -pkin(7) * t132 + t35;
t27 = pkin(7) * t96 + t33;
t18 = t20 * t117;
t5 = -t122 * t23 + t172 * t22;
t4 = -qJD(5) * t12 - t122 * t29 + t172 * t27;
t3 = qJD(5) * t11 + t122 * t27 + t172 * t29;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, -t156, 0, t97 * t123 - t98 * t124 + (-t104 * t123 + t105 * t124) * qJD(3), 0, 0, 0, 0, 0, 0, -t82, t155, -t130 + t164, t25 * t100 - t132 * t30 - t24 * t99 - t31 * t96, 0, 0, 0, 0, 0, 0, -t18, t161, t135 + t175, -t128 * t44 - t19 * t6 - t2 * t43 - t20 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t136, t154 * t180, t156, -0.2e1 * t136, -t157, 0, -pkin(6) * t156 + t123 * t137, pkin(6) * t157 + t124 * t137, t127, t127 * pkin(6), t100 * t81 - t131 * t96, t130 + t164, -t155, t132 * t92 + t80 * t99, -t82, 0, t115 * t80 + (t106 * t100 + t33 + (qJD(2) * t99 + t92) * t170) * qJD(3), -t106 * t96 + t115 * t81 + (t170 * t179 - t35) * qJD(3), -t24 * t100 - t131 * t33 - t132 * t31 - t25 * t99 + t30 * t96 - t35 * t92 - t54 * t81 - t55 * t80, t24 * t54 + t25 * t55 + t30 * t33 + t31 * t35 + (t106 + t153) * t123 * t163, t133 * t44 + t134 * t19, -t135 + t175, -t161, -t129 * t43 - t20 * t41, -t18, 0, t117 * t4 - t129 * t64 + t20 * t50 - t41 * t57 + t43 * t51, -t117 * t3 + t133 * t64 - t134 * t57 - t19 * t50 + t44 * t51, -t11 * t133 + t12 * t129 + t128 * t43 + t134 * t4 + t19 * t5 - t2 * t44 - t20 * t6 + t3 * t41, t11 * t2 - t12 * t128 + t3 * t6 + t4 * t5 + t50 * t57 + t51 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t154 * t126, 0, t144, 0, 0, t126 * pkin(2) * t123, pkin(2) * t158 - t114 + (t104 + t145) * qJD(3), 0, 0, t166, -t90 + t178, -t109 + (t92 + t143) * qJD(3), -t166, 0, 0, -qJD(3) * t32 - t106 * t131 - t146 * t92 + t24, qJD(3) * t34 + t106 * t92 - t131 * t146 - t25, (t31 + t32) * t131 + (-t30 + t34) * t92 + (-t120 * t80 - t121 * t81) * pkin(3), -t30 * t32 - t31 * t34 + (-t106 * t152 + t120 * t25 + t121 * t24) * pkin(3), t167, t183, t184, -t167, t185, 0, -t117 * t174 + t41 * t56 + t182, -t117 * t173 + t134 * t56 + t181, -t133 * t84 + t85 * t129 + (t173 + t5) * t41 + (-t174 - t6) * t134, -t128 * t85 + t173 * t6 - t174 * t5 + t2 * t84 - t50 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179 * qJD(3), -t109 + (-t92 + t143) * qJD(3), -t90 - t178, t131 * t30 + t31 * t92 + t111, 0, 0, 0, 0, 0, 0, -t129 - t160, t133 + t159, -t168 - t169, -t134 * t5 - t41 * t6 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, t183, t184, -t167, t185, 0, t6 * t117 + t182, t117 * t5 + t181, 0, 0;];
tauc_reg = t1;
