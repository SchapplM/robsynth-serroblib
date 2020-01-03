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
% tauc_reg [4x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 17:04:54
% EndTime: 2019-12-31 17:04:56
% DurationCPUTime: 0.57s
% Computational Cost: add. (712->129), mult. (1987->200), div. (0->0), fcn. (1431->6), ass. (0->96)
t88 = sin(qJ(4));
t108 = qJD(4) * t88;
t87 = cos(pkin(7));
t91 = cos(qJ(2));
t115 = t87 * t91;
t86 = sin(pkin(7));
t89 = sin(qJ(2));
t70 = -t86 * t89 + t115;
t60 = t70 * qJD(1);
t90 = cos(qJ(4));
t51 = t90 * t60;
t71 = t86 * t91 + t87 * t89;
t61 = t71 * qJD(2);
t54 = qJD(1) * t61;
t106 = qJD(1) * qJD(2);
t102 = t89 * t106;
t55 = -t86 * t102 + t106 * t115;
t62 = t71 * qJD(1);
t1 = qJD(4) * t51 - t62 * t108 - t88 * t54 + t90 * t55;
t24 = -t88 * t62 + t51;
t83 = qJD(2) + qJD(4);
t117 = t24 * t83;
t130 = t1 - t117;
t96 = t88 * t60 + t90 * t62;
t129 = t96 * t24;
t118 = t96 * t83;
t2 = t96 * qJD(4) + t90 * t54 + t88 * t55;
t128 = -t2 + t118;
t127 = -t24 ^ 2 + t96 ^ 2;
t120 = t60 * pkin(6);
t111 = -qJ(3) - pkin(5);
t79 = t111 * t91;
t76 = qJD(1) * t79;
t116 = t87 * t76;
t109 = qJD(2) * pkin(2);
t78 = t111 * t89;
t75 = qJD(1) * t78;
t69 = t75 + t109;
t27 = t86 * t69 - t116;
t12 = t27 + t120;
t104 = -t91 * pkin(2) - pkin(1);
t98 = t104 * qJD(1);
t77 = qJD(3) + t98;
t32 = -t60 * pkin(3) + t77;
t126 = t12 * t108 - t32 * t24;
t124 = -0.2e1 * t106;
t123 = qJD(4) - t83;
t100 = qJD(2) * t111;
t58 = t91 * qJD(3) + t89 * t100;
t44 = t58 * qJD(1);
t59 = -t89 * qJD(3) + t91 * t100;
t45 = t59 * qJD(1);
t13 = -t86 * t44 + t87 * t45;
t6 = -t55 * pkin(6) + t13;
t14 = t87 * t44 + t86 * t45;
t7 = -t54 * pkin(6) + t14;
t122 = -t32 * t96 + t90 * t6 - t88 * t7;
t121 = pkin(2) * t86;
t119 = t62 * pkin(6);
t65 = t86 * t76;
t93 = qJD(1) ^ 2;
t114 = t91 * t93;
t92 = qJD(2) ^ 2;
t113 = t92 * t89;
t112 = t92 * t91;
t18 = t87 * t58 + t86 * t59;
t31 = t87 * t75 + t65;
t35 = t86 * t78 - t87 * t79;
t110 = t89 ^ 2 - t91 ^ 2;
t107 = t89 * qJD(1);
t105 = t89 * t109;
t17 = -t86 * t58 + t87 * t59;
t26 = t87 * t69 + t65;
t30 = -t86 * t75 + t116;
t34 = t87 * t78 + t86 * t79;
t99 = pkin(1) * t124;
t11 = qJD(2) * pkin(3) - t119 + t26;
t97 = -t88 * t11 - t90 * t12;
t95 = t90 * t70 - t88 * t71;
t29 = t88 * t70 + t90 * t71;
t82 = t87 * pkin(2) + pkin(3);
t81 = pkin(2) * t102;
t64 = t70 * qJD(2);
t47 = -t70 * pkin(3) + t104;
t40 = t61 * pkin(3) + t105;
t39 = pkin(2) * t107 + t62 * pkin(3);
t33 = t54 * pkin(3) + t81;
t20 = t70 * pkin(6) + t35;
t19 = -t71 * pkin(6) + t34;
t16 = t31 - t119;
t15 = t30 - t120;
t9 = -t61 * pkin(6) + t18;
t8 = -t64 * pkin(6) + t17;
t4 = t29 * qJD(4) + t90 * t61 + t88 * t64;
t3 = t95 * qJD(4) - t88 * t61 + t90 * t64;
t5 = [0, 0, 0, 0.2e1 * t91 * t102, t110 * t124, t112, -t113, 0, -pkin(5) * t112 + t89 * t99, pkin(5) * t113 + t91 * t99, -t13 * t71 + t14 * t70 - t17 * t62 + t18 * t60 - t26 * t64 - t27 * t61 - t34 * t55 - t35 * t54, t13 * t34 + t14 * t35 + t26 * t17 + t27 * t18 + (t77 + t98) * t105, t1 * t29 + t3 * t96, t1 * t95 - t29 * t2 + t24 * t3 - t4 * t96, t3 * t83, -t4 * t83, 0, -t40 * t24 + t47 * t2 - t33 * t95 + t32 * t4 + (t90 * t8 - t88 * t9 + (-t19 * t88 - t20 * t90) * qJD(4)) * t83, t40 * t96 + t47 * t1 + t33 * t29 + t32 * t3 - (t88 * t8 + t90 * t9 + (t19 * t90 - t20 * t88) * qJD(4)) * t83; 0, 0, 0, -t89 * t114, t110 * t93, 0, 0, 0, t93 * pkin(1) * t89, pkin(1) * t114, (t27 + t30) * t62 + (t26 - t31) * t60 + (-t54 * t86 - t55 * t87) * pkin(2), -t26 * t30 - t27 * t31 + (-t77 * t107 + t13 * t87 + t14 * t86) * pkin(2), -t129, t127, t130, t128, 0, t39 * t24 - (t90 * t15 - t88 * t16) * t83 + ((-t90 * t121 - t82 * t88) * t83 + t97) * qJD(4) + t122, -t90 * t7 - t88 * t6 - t39 * t96 + (t88 * t15 + t90 * t16) * t83 + (-(-t88 * t121 + t82 * t90) * t83 - t90 * t11) * qJD(4) + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60 ^ 2 - t62 ^ 2, t26 * t62 - t27 * t60 + t81, 0, 0, 0, 0, 0, t2 + t118, t1 + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t127, t130, t128, 0, t123 * t97 + t122, (-t12 * t83 - t6) * t88 + (-t123 * t11 - t7) * t90 + t126;];
tauc_reg = t5;
