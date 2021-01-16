% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:22
% EndTime: 2021-01-15 10:46:27
% DurationCPUTime: 0.71s
% Computational Cost: add. (780->143), mult. (2181->221), div. (0->0), fcn. (1551->6), ass. (0->102)
t92 = sin(qJ(4));
t112 = qJD(4) * t92;
t90 = sin(pkin(7));
t91 = cos(pkin(7));
t93 = sin(qJ(2));
t95 = cos(qJ(2));
t72 = t90 * t95 + t91 * t93;
t115 = qJD(1) * t72;
t122 = t91 * t95;
t108 = qJD(1) * t122;
t113 = qJD(1) * t93;
t60 = t90 * t113 - t108;
t94 = cos(qJ(4));
t51 = t94 * t60;
t62 = t72 * qJD(2);
t54 = qJD(1) * t62;
t111 = qJD(1) * qJD(2);
t105 = t95 * t111;
t106 = t93 * t111;
t81 = t90 * t106;
t55 = t91 * t105 - t81;
t1 = -qJD(4) * t51 - t112 * t115 - t92 * t54 + t94 * t55;
t24 = -t115 * t92 - t51;
t87 = qJD(2) + qJD(4);
t124 = t24 * t87;
t137 = t1 - t124;
t99 = -t115 * t94 + t92 * t60;
t125 = t99 * t87;
t98 = t99 * qJD(4) - t94 * t54 - t92 * t55;
t136 = t98 - t125;
t135 = t99 * t24;
t134 = -t24 ^ 2 + t99 ^ 2;
t127 = t60 * pkin(6);
t118 = -qJ(3) - pkin(5);
t80 = t118 * t95;
t77 = qJD(1) * t80;
t123 = t91 * t77;
t116 = qJD(2) * pkin(2);
t79 = t118 * t93;
t76 = qJD(1) * t79;
t70 = t76 + t116;
t27 = t90 * t70 - t123;
t12 = t27 - t127;
t86 = -t95 * pkin(2) - pkin(1);
t114 = qJD(1) * t86;
t78 = qJD(3) + t114;
t32 = t60 * pkin(3) + t78;
t133 = t12 * t112 - t32 * t24;
t103 = qJD(2) * t118;
t58 = t95 * qJD(3) + t93 * t103;
t44 = t58 * qJD(1);
t59 = -t93 * qJD(3) + t95 * t103;
t45 = t59 * qJD(1);
t13 = -t90 * t44 + t91 * t45;
t6 = -t55 * pkin(6) + t13;
t14 = t91 * t44 + t90 * t45;
t7 = -t54 * pkin(6) + t14;
t132 = t32 * t99 + t94 * t6 - t92 * t7;
t131 = -0.2e1 * t111;
t130 = qJD(4) - t87;
t129 = pkin(2) * t90;
t128 = pkin(2) * t93;
t126 = t115 * pkin(6);
t66 = t90 * t77;
t97 = qJD(1) ^ 2;
t121 = t95 * t97;
t96 = qJD(2) ^ 2;
t120 = t96 * t93;
t119 = t96 * t95;
t18 = t91 * t58 + t90 * t59;
t31 = t91 * t76 + t66;
t35 = t90 * t79 - t91 * t80;
t117 = t93 ^ 2 - t95 ^ 2;
t110 = t93 * t116;
t109 = pkin(2) * t113;
t17 = -t90 * t58 + t91 * t59;
t26 = t91 * t70 + t66;
t30 = -t90 * t76 + t123;
t34 = t91 * t79 + t90 * t80;
t102 = 0.2e1 * t115;
t101 = pkin(1) * t131;
t11 = qJD(2) * pkin(3) - t126 + t26;
t100 = -t92 * t11 - t94 * t12;
t71 = t90 * t93 - t122;
t28 = t94 * t71 + t92 * t72;
t29 = -t92 * t71 + t94 * t72;
t85 = t91 * pkin(2) + pkin(3);
t83 = pkin(2) * t106;
t65 = t71 * qJD(2);
t47 = t71 * pkin(3) + t86;
t40 = t62 * pkin(3) + t110;
t39 = pkin(3) * t115 + t109;
t33 = t54 * pkin(3) + t83;
t20 = -t71 * pkin(6) + t35;
t19 = -t72 * pkin(6) + t34;
t16 = t31 - t126;
t15 = t30 + t127;
t9 = -t62 * pkin(6) + t18;
t8 = t65 * pkin(6) + t17;
t4 = t29 * qJD(4) + t94 * t62 - t92 * t65;
t3 = -t28 * qJD(4) - t92 * t62 - t94 * t65;
t2 = [0, 0, 0, 0.2e1 * t93 * t105, t117 * t131, t119, -t120, 0, -pkin(5) * t119 + t101 * t93, pkin(5) * t120 + t95 * t101, t86 * t54 + t78 * t62 + (t17 + (qJD(1) * t71 + t60) * t128) * qJD(2), t86 * t55 - t78 * t65 + (t102 * t128 - t18) * qJD(2), -t115 * t17 - t13 * t72 - t14 * t71 - t18 * t60 + t26 * t65 - t27 * t62 - t34 * t55 - t35 * t54, t13 * t34 + t14 * t35 + t26 * t17 + t27 * t18 + (t78 + t114) * t110, t1 * t29 - t3 * t99, -t1 * t28 + t24 * t3 + t29 * t98 + t4 * t99, t3 * t87, -t4 * t87, 0, -t40 * t24 - t47 * t98 + t33 * t28 + t32 * t4 + (t94 * t8 - t92 * t9 + (-t19 * t92 - t20 * t94) * qJD(4)) * t87, -t40 * t99 + t47 * t1 + t33 * t29 + t32 * t3 - (t92 * t8 + t94 * t9 + (t19 * t94 - t20 * t92) * qJD(4)) * t87; 0, 0, 0, -t93 * t121, t117 * t97, 0, 0, 0, t97 * pkin(1) * t93, pkin(1) * t121, -t30 * qJD(2) - t109 * t60 - t115 * t78 + t13, t31 * qJD(2) - t109 * t115 + t78 * t60 - t14, (t27 + t30) * t115 + (-t26 + t31) * t60 + (-t54 * t90 - t55 * t91) * pkin(2), -t26 * t30 - t27 * t31 + (-t113 * t78 + t13 * t91 + t14 * t90) * pkin(2), t135, t134, t137, t136, 0, t39 * t24 - (t94 * t15 - t92 * t16) * t87 + ((-t94 * t129 - t85 * t92) * t87 + t100) * qJD(4) + t132, -t94 * t7 - t92 * t6 + t39 * t99 + (t92 * t15 + t94 * t16) * t87 + (-(-t92 * t129 + t85 * t94) * t87 - t94 * t11) * qJD(4) + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * qJD(2), -t81 + (-t60 + t108) * qJD(2), -t115 ^ 2 - t60 ^ 2, t115 * t26 + t27 * t60 + t83, 0, 0, 0, 0, 0, -t98 - t125, t1 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t134, t137, t136, 0, t130 * t100 + t132, (-t12 * t87 - t6) * t92 + (-t130 * t11 - t7) * t94 + t133;];
tauc_reg = t2;
