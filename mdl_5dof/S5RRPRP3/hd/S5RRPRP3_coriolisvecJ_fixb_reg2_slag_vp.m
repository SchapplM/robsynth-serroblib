% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:13
% EndTime: 2019-12-31 19:51:17
% DurationCPUTime: 1.02s
% Computational Cost: add. (2122->186), mult. (3713->221), div. (0->0), fcn. (2454->6), ass. (0->118)
t103 = cos(qJ(2));
t147 = pkin(1) * qJD(1);
t116 = -t103 * t147 + qJD(3);
t100 = cos(pkin(8));
t99 = sin(pkin(8));
t148 = t100 ^ 2 + t99 ^ 2;
t146 = pkin(1) * qJD(2);
t135 = qJD(1) * t146;
t98 = qJD(1) + qJD(2);
t77 = t98 * qJD(3) + t103 * t135;
t173 = t148 * t77;
t161 = cos(qJ(4));
t130 = t161 * t100;
t101 = sin(qJ(4));
t144 = t101 * t99;
t110 = t130 - t144;
t84 = (-pkin(7) - qJ(3)) * t99;
t94 = t100 * pkin(7);
t85 = t100 * qJ(3) + t94;
t111 = -t101 * t85 + t161 * t84;
t154 = t111 * qJD(4) + t116 * t110;
t79 = t101 * t100 + t161 * t99;
t170 = t79 * t98;
t171 = t79 * t77;
t93 = -t100 * pkin(3) - pkin(2);
t58 = t93 * t98 + t116;
t139 = t98 * t144;
t59 = -t98 * t130 + t139;
t26 = t59 * pkin(4) - qJ(5) * t170 + t58;
t172 = t26 * t170 + t171;
t165 = t170 ^ 2;
t55 = t59 ^ 2;
t169 = -t55 - t165;
t168 = -t55 + t165;
t102 = sin(qJ(2));
t167 = t100 * t102;
t166 = t154 * qJD(4);
t92 = t102 * pkin(1) + qJ(3);
t74 = (-pkin(7) - t92) * t99;
t75 = t100 * t92 + t94;
t112 = -t101 * t75 + t161 * t74;
t137 = t102 * t147;
t81 = t98 * qJ(3) + t137;
t133 = pkin(7) * t98 + t81;
t49 = t133 * t99;
t50 = t133 * t100;
t32 = -t101 * t49 + t161 * t50;
t9 = t32 * qJD(4) + t171;
t164 = t9 * t112;
t163 = t9 * t111;
t162 = t9 * t79;
t160 = t103 * pkin(1);
t158 = t170 * t59;
t134 = qJD(4) * t144;
t127 = qJD(4) * t161;
t120 = t100 * t127;
t82 = t98 * t120;
t53 = t98 * t134 - t82;
t73 = t79 * qJD(4);
t54 = t98 * t73;
t90 = t102 * t135;
t113 = t54 * pkin(4) + t53 * qJ(5) + t90;
t15 = -qJD(5) * t170 + t113;
t157 = -t110 * t15 + t26 * t73;
t72 = -t120 + t134;
t156 = -t15 * t79 + t26 * t72;
t29 = qJD(4) * qJ(5) + t32;
t155 = t29 - t32;
t48 = t101 * t84 + t161 * t85;
t153 = t48 * qJD(4) + t116 * t79;
t152 = -t110 * t90 + t58 * t73;
t151 = -t49 * t127 + t77 * t130;
t150 = -t58 * t72 + t79 * t90;
t145 = t101 * t50;
t143 = qJD(4) * t73;
t88 = t103 * t146 + qJD(3);
t23 = t112 * qJD(4) + t110 * t88;
t142 = t23 * qJD(4);
t44 = t101 * t74 + t161 * t75;
t24 = t44 * qJD(4) + t79 * t88;
t141 = t24 * qJD(4);
t31 = -t161 * t49 - t145;
t140 = qJD(5) - t31;
t138 = t102 * t146;
t131 = t148 * t88;
t8 = (-qJD(4) * t50 - t77 * t99) * t101 + t151;
t129 = t110 * t8 + t31 * t72 - t32 * t73 + t162;
t28 = -qJD(4) * pkin(4) + t140;
t114 = t77 * t144 - t151;
t7 = (qJD(5) - t145) * qJD(4) - t114;
t128 = t110 * t7 - t28 * t72 - t29 * t73 + t162;
t126 = t148 * t103;
t125 = t31 + t145;
t124 = t153 * qJD(4);
t123 = t148 * qJD(3);
t122 = t98 * t138;
t121 = t98 * t137;
t33 = t73 * pkin(4) + t72 * qJ(5) - t79 * qJD(5);
t119 = -t33 + t137;
t118 = (-qJD(2) + t98) * t147;
t117 = (-qJD(1) - t98) * t146;
t115 = -t110 * t54 + t59 * t73;
t109 = t112 * t53 + t170 * t24 - t23 * t59 - t44 * t54;
t45 = -pkin(4) * t110 - t79 * qJ(5) + t93;
t2 = -t110 * t53 - t170 * t73 - t79 * t54 + t72 * t59;
t106 = t111 * t53 + t153 * t170 - t154 * t59 - t48 * t54;
t105 = 0.2e1 * t170 * qJD(4);
t86 = t99 * t90;
t83 = t93 - t160;
t80 = -t98 * pkin(2) + t116;
t67 = t72 * qJD(4);
t40 = t45 - t160;
t39 = pkin(4) * t170 + t59 * qJ(5);
t38 = t82 + (t59 - t139) * qJD(4);
t37 = -t82 + (t59 + t139) * qJD(4);
t30 = t33 + t138;
t18 = -t170 * t72 - t53 * t79;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90 - t122, t103 * t117, 0, 0, 0, 0, 0, 0, 0, 0, t117 * t167, t99 * t122 + t86, t98 * t131 + t173, t81 * t131 + t92 * t173 + (t80 + (-pkin(2) - t160) * qJD(1)) * t138, t18, t2, -t67, t115, -t143, 0, t59 * t138 + t83 * t54 - t141 + t152, t138 * t170 - t83 * t53 - t142 + t150, t109 + t129, t32 * t23 - t31 * t24 - t164 + t8 * t44 + (qJD(1) * t83 + t58) * t138, t18, -t67, -t2, 0, t143, t115, t30 * t59 + t40 * t54 - t141 + t157, t109 + t128, -t170 * t30 + t40 * t53 + t142 + t156, t15 * t40 + t29 * t23 + t28 * t24 + t26 * t30 + t7 * t44 - t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90 + t121, t103 * t118, 0, 0, 0, 0, 0, 0, 0, 0, t118 * t167, -t99 * t121 + t86, (-t126 * t147 + t123) * t98 + t173, t81 * t123 + qJ(3) * t173 + (-t81 * t126 + (-pkin(2) * qJD(2) - t80) * t102) * t147, t18, t2, -t67, t115, -t143, 0, -t59 * t137 + t93 * t54 - t124 + t152, -t137 * t170 - t93 * t53 + t150 - t166, t106 + t129, -t163 + t8 * t48 + t154 * t32 - t153 * t31 + (qJD(2) * t93 - t58) * t137, t18, -t67, -t2, 0, t143, t115, -t119 * t59 + t45 * t54 - t124 + t157, t106 + t128, t119 * t170 + t45 * t53 + t156 + t166, -t119 * t26 + t15 * t45 + t153 * t28 + t154 * t29 + t7 * t48 - t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148 * t98 ^ 2, -t148 * t98 * t81 + t90, 0, 0, 0, 0, 0, 0, t105, -t37, t169, t170 * t31 + t32 * t59 + t90, 0, 0, 0, 0, 0, 0, t105, t169, t37, t29 * t59 + (-qJD(5) - t28) * t170 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t168, t38, -t158, 0, 0, -t170 * t58 - t171, t125 * qJD(4) + t58 * t59 + t114, 0, 0, t158, t38, -t168, 0, 0, -t158, -t39 * t59 - t172, pkin(4) * t53 - t54 * qJ(5) + t155 * t170 + (t28 - t140) * t59, -t26 * t59 + t39 * t170 + (0.2e1 * qJD(5) - t125) * qJD(4) - t114, -t9 * pkin(4) + t7 * qJ(5) + t140 * t29 - t26 * t39 - t28 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t38, -qJD(4) ^ 2 - t165, -t155 * qJD(4) + t172;];
tauc_reg = t1;
