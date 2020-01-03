% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [4x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:14
% EndTime: 2019-12-31 17:26:18
% DurationCPUTime: 1.12s
% Computational Cost: add. (1231->168), mult. (3180->261), div. (0->0), fcn. (2234->6), ass. (0->109)
t131 = cos(qJ(3));
t100 = qJD(1) * t131;
t76 = sin(qJ(2));
t112 = qJD(1) * t76;
t75 = sin(qJ(3));
t78 = cos(qJ(2));
t137 = t78 * t100 - t75 * t112;
t71 = qJD(2) + qJD(3);
t27 = t137 * t71;
t108 = qJD(1) * qJD(2);
t136 = -0.2e1 * t108;
t119 = t75 * t78;
t52 = t131 * t76 + t119;
t135 = qJD(1) * t52;
t43 = qJD(4) - t137;
t134 = -qJD(4) + t43;
t74 = sin(qJ(4));
t47 = -qJD(1) * t119 - t76 * t100;
t77 = cos(qJ(4));
t89 = t77 * t47 - t74 * t71;
t14 = -t89 * qJD(4) + t74 * t27;
t111 = qJD(3) * t75;
t113 = qJD(2) * pkin(2);
t132 = pkin(5) + pkin(6);
t58 = t132 * t76;
t53 = qJD(1) * t58;
t48 = -t53 + t113;
t103 = qJD(2) * t132;
t91 = qJD(1) * t103;
t49 = t76 * t91;
t50 = t78 * t91;
t59 = t132 * t78;
t55 = qJD(1) * t59;
t99 = t131 * qJD(3);
t10 = -t55 * t111 - t131 * t49 + t48 * t99 - t75 * t50;
t101 = t131 * t50 - t75 * t49;
t104 = t131 * t55;
t32 = t75 * t48 + t104;
t11 = t32 * qJD(3) + t101;
t54 = t76 * t103;
t56 = t78 * t103;
t84 = -t131 * t58 - t75 * t59;
t15 = t84 * qJD(3) - t131 * t54 - t75 * t56;
t70 = -t78 * pkin(2) - pkin(1);
t57 = t70 * qJD(1);
t22 = -pkin(3) * t137 + t47 * pkin(7) + t57;
t120 = t75 * t55;
t31 = t131 * t48 - t120;
t24 = -t71 * pkin(3) - t31;
t36 = t71 * t52;
t28 = t36 * qJD(1);
t83 = t131 * t78 - t75 * t76;
t30 = -pkin(3) * t83 - t52 * pkin(7) + t70;
t35 = t71 * t83;
t41 = t131 * t59 - t75 * t58;
t133 = t11 * t52 + t24 * t35 - t41 * t28 - (qJD(4) * t30 + t15) * t43 + (qJD(4) * t22 + t10) * t83;
t109 = qJD(4) * t77;
t110 = qJD(4) * t74;
t13 = t71 * t109 + t47 * t110 + t77 * t27;
t130 = t13 * t74;
t129 = t24 * t137;
t128 = t24 * t52;
t127 = t30 * t28;
t37 = -t74 * t47 - t77 * t71;
t126 = t37 * t43;
t125 = t89 * t43;
t124 = t43 * t47;
t123 = t47 * t137;
t121 = t74 * t28;
t118 = t77 * t28;
t80 = qJD(1) ^ 2;
t117 = t78 * t80;
t79 = qJD(2) ^ 2;
t116 = t79 * t76;
t115 = t79 * t78;
t114 = t76 ^ 2 - t78 ^ 2;
t107 = t76 * t113;
t106 = pkin(2) * t112;
t102 = t76 * t108;
t98 = t43 * t77;
t29 = -t47 * pkin(3) - pkin(7) * t137;
t68 = t75 * pkin(2) + pkin(7);
t95 = qJD(4) * t68 + t106 + t29;
t94 = pkin(1) * t136;
t25 = t71 * pkin(7) + t32;
t6 = t74 * t22 + t77 * t25;
t93 = t24 * t109 + t11 * t74 - t6 * t47;
t33 = -t75 * t53 + t104;
t92 = pkin(2) * t111 - t33;
t5 = t77 * t22 - t74 * t25;
t90 = -t68 * t28 - t129;
t88 = -t11 * t77 + t24 * t110 + t5 * t47;
t34 = -t131 * t53 - t120;
t87 = -pkin(2) * t99 + t34;
t86 = t57 * t47 - t101;
t85 = -t52 * t110 + t77 * t35;
t81 = -t137 * t57 - t10;
t69 = -t131 * pkin(2) - pkin(3);
t19 = -t137 ^ 2 + t47 ^ 2;
t18 = (-t135 - t47) * t71;
t16 = t41 * qJD(3) + t131 * t56 - t75 * t54;
t12 = t36 * pkin(3) - t35 * pkin(7) + t107;
t8 = pkin(2) * t102 + t28 * pkin(3) - t27 * pkin(7);
t7 = t77 * t8;
t4 = t43 * t98 - t47 * t89 + t121;
t3 = -t43 ^ 2 * t74 - t37 * t47 + t118;
t2 = -t89 * t98 + t130;
t1 = (t13 - t126) * t77 + (-t14 + t125) * t74;
t9 = [0, 0, 0, 0.2e1 * t78 * t102, t114 * t136, t115, -t116, 0, -pkin(5) * t115 + t76 * t94, pkin(5) * t116 + t78 * t94, t27 * t52 - t47 * t35, t137 * t35 + t27 * t83 - t52 * t28 + t47 * t36, t35 * t71, -t36 * t71, 0, -t16 * t71 + t70 * t28 + t57 * t36 + (-qJD(1) * t83 - t137) * t107, -t15 * t71 + t70 * t27 + t57 * t35 + (-t47 + t135) * t107, t13 * t77 * t52 - t85 * t89, (-t37 * t77 + t74 * t89) * t35 + (-t130 - t14 * t77 + (t37 * t74 + t77 * t89) * qJD(4)) * t52, t52 * t118 - t13 * t83 - t36 * t89 + t85 * t43, -t52 * t121 + t14 * t83 - t37 * t36 + (-t52 * t109 - t74 * t35) * t43, -t28 * t83 + t43 * t36, -t84 * t14 + t16 * t37 + t5 * t36 - t7 * t83 + (t12 * t43 + t127 + (t25 * t83 - t41 * t43 + t128) * qJD(4)) * t77 + t133 * t74, -t84 * t13 - t16 * t89 - t6 * t36 + (-(-qJD(4) * t41 + t12) * t43 - t127 + (-qJD(4) * t25 + t8) * t83 - qJD(4) * t128) * t74 + t133 * t77; 0, 0, 0, -t76 * t117, t114 * t80, 0, 0, 0, t80 * pkin(1) * t76, pkin(1) * t117, t123, t19, 0, t18, 0, t137 * t106 + t33 * t71 + (-t104 + (-pkin(2) * t71 - t48) * t75) * qJD(3) + t86, t34 * t71 + (t47 * t112 - t71 * t99) * pkin(2) + t81, t2, t1, t4, t3, t124, t69 * t14 + t90 * t74 + t92 * t37 + (t87 * t74 - t95 * t77) * t43 + t88, t69 * t13 + t90 * t77 - t92 * t89 + (t95 * t74 + t87 * t77) * t43 + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t19, 0, t18, 0, t86 + (-qJD(3) + t71) * t32, t31 * t71 + t81, t2, t1, t4, t3, t124, -pkin(3) * t14 - (t77 * t29 - t74 * t31) * t43 - t32 * t37 - t74 * t129 + (-t43 * t109 - t121) * pkin(7) + t88, -pkin(3) * t13 + (t74 * t29 + t77 * t31) * t43 + t32 * t89 - t77 * t129 + (t43 * t110 - t118) * pkin(7) + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89 * t37, -t37 ^ 2 + t89 ^ 2, t13 + t126, -t14 - t125, t28, -t74 * t10 + t134 * t6 + t24 * t89 + t7, -t77 * t10 + t134 * t5 + t24 * t37 - t74 * t8;];
tauc_reg = t9;
