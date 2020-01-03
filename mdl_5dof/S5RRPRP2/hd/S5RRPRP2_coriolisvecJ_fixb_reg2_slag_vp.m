% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP2
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:43
% EndTime: 2019-12-31 19:49:45
% DurationCPUTime: 0.62s
% Computational Cost: add. (1130->144), mult. (2360->179), div. (0->0), fcn. (1290->6), ass. (0->107)
t107 = pkin(1) * qJD(1);
t68 = sin(qJ(2));
t100 = t68 * t107;
t65 = sin(pkin(8));
t51 = t65 * t100;
t66 = cos(pkin(8));
t70 = cos(qJ(2));
t99 = t70 * t107;
t37 = t66 * t99 - t51;
t67 = sin(qJ(4));
t63 = t67 ^ 2;
t69 = cos(qJ(4));
t64 = t69 ^ 2;
t132 = t63 + t64;
t102 = t69 * qJD(3);
t62 = qJD(1) + qJD(2);
t47 = t62 * pkin(2) + t99;
t25 = t66 * t100 + t65 * t47;
t19 = t62 * pkin(7) + t25;
t116 = t67 * t19;
t13 = t102 - t116;
t131 = qJD(5) - t13;
t30 = t37 * qJD(2);
t111 = qJD(4) * t102 + t69 * t30;
t2 = (qJD(5) - t116) * qJD(4) + t111;
t115 = t67 * t30;
t14 = t67 * qJD(3) + t69 * t19;
t6 = t14 * qJD(4) + t115;
t4 = t6 * t67;
t130 = t2 * t69 + t4;
t105 = qJD(4) * t67;
t5 = -t19 * t105 + t111;
t129 = t5 * t69 + t4 + (-t13 * t69 - t14 * t67) * qJD(4);
t11 = -qJD(4) * pkin(4) + t131;
t101 = qJD(4) * qJ(5);
t12 = t14 + t101;
t108 = -t63 + t64;
t128 = 0.2e1 * t108 * t62 * qJD(4);
t127 = t6 * t69;
t126 = t66 * pkin(2);
t117 = t66 * t68;
t57 = t70 * pkin(1) + pkin(2);
t109 = pkin(1) * t117 + t65 * t57;
t34 = pkin(7) + t109;
t71 = qJD(4) ^ 2;
t125 = t34 * t71;
t83 = t65 * t70 + t117;
t78 = pkin(1) * t83;
t35 = qJD(1) * t78;
t124 = t35 * t62;
t36 = qJD(2) * t78;
t123 = t36 * t62;
t122 = t37 * t62;
t106 = pkin(1) * qJD(2);
t118 = t65 * t68;
t38 = (t66 * t70 - t118) * t106;
t121 = t38 * t62;
t54 = t65 * pkin(2) + pkin(7);
t120 = t54 * t71;
t119 = t62 * t69;
t104 = qJD(4) * t69;
t24 = t66 * t47 - t51;
t18 = -t62 * pkin(3) - t24;
t29 = qJD(1) * t36;
t114 = t18 * t104 + t29 * t67;
t113 = t132 * t121;
t112 = t37 * t105 + t35 * t119;
t103 = t67 * qJD(5);
t40 = pkin(4) * t105 - t69 * t101 - t103;
t110 = t35 - t40;
t98 = t62 * t105;
t87 = pkin(4) * t67 - qJ(5) * t69;
t8 = t29 + (t87 * qJD(4) - t103) * t62;
t97 = -t8 - t120;
t96 = t13 + t116;
t81 = -t69 * pkin(4) - t67 * qJ(5) - pkin(3);
t94 = -pkin(1) * t118 + t66 * t57;
t20 = t81 - t94;
t95 = t20 * t62 - t38;
t91 = t69 * t98;
t90 = t132 * t122;
t89 = (-qJD(2) + t62) * t107;
t88 = (-qJD(1) - t62) * t106;
t86 = t11 * t67 + t12 * t69;
t85 = t13 * t67 - t14 * t69;
t84 = t123 + t125;
t33 = -pkin(3) - t94;
t82 = qJD(4) * (t33 * t62 - t38);
t15 = t36 + t40;
t80 = -t15 * t62 - t125 - t8;
t79 = t11 * t104 - t12 * t105 + t130;
t76 = t83 * qJD(1) * t106;
t73 = (t11 * t69 - t12 * t67) * qJD(4) + t130;
t61 = t62 ^ 2;
t60 = t71 * t69;
t59 = t71 * t67;
t55 = -pkin(3) - t126;
t49 = t67 * t61 * t69;
t46 = -0.2e1 * t91;
t45 = 0.2e1 * t91;
t44 = t108 * t61;
t43 = t81 - t126;
t41 = t87 * t62;
t16 = t18 * t105;
t10 = t81 * t62 - t24;
t7 = t10 * t105;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t88, t70 * t88, 0, 0, 0, 0, 0, 0, 0, 0, -t76 - t123, -t30 - t121, 0, t30 * t109 - t24 * t36 + t25 * t38 - t29 * t94, t45, t128, t60, t46, -t59, 0, t16 + t67 * t82 + (-t29 - t84) * t69, t84 * t67 + t69 * t82 + t114, t129 + t113, t129 * t34 + t18 * t36 + t29 * t33 - t85 * t38, t45, t60, -t128, 0, t59, t46, t95 * t105 + t80 * t69 + t7, t79 + t113, t80 * t67 + (-t10 - t95) * t104, t10 * t15 + t8 * t20 + t73 * t34 + t86 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t89, t70 * t89, 0, 0, 0, 0, 0, 0, 0, 0, -t76 + t124, -t30 + t122, 0, t24 * t35 - t25 * t37 + (-t29 * t66 + t30 * t65) * pkin(2), t45, t128, t60, t46, -t59, 0, t55 * t98 + t16 + (-t29 - t120) * t69 + t112, (t120 - t124) * t67 + (t55 * t62 + t37) * t104 + t114, -t90 + t129, t129 * t54 - t18 * t35 + t29 * t55 + t85 * t37, t45, t60, -t128, 0, t59, t46, t43 * t98 + t7 + (-t40 * t62 + t97) * t69 + t112, -t90 + t79, (-t43 * t62 - t10 - t37) * t104 + (t110 * t62 + t97) * t67, -t110 * t10 - t86 * t37 + t8 * t43 + t73 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t60, 0, -t85 * qJD(4) + t5 * t67 - t127, 0, 0, 0, 0, 0, 0, -t59, 0, t60, t86 * qJD(4) + t2 * t67 - t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t44, 0, t49, 0, 0, (-t18 * t62 - t30) * t67, t96 * qJD(4) - t18 * t119 - t111, 0, 0, -t49, 0, t44, 0, 0, t49, -t115 + (-t10 * t67 + t41 * t69) * t62, 0, (t10 * t69 + t41 * t67) * t62 + (0.2e1 * qJD(5) - t96) * qJD(4) + t111, -t6 * pkin(4) + t2 * qJ(5) - t10 * t41 - t11 * t14 + t131 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, -t63 * t61 - t71, (t10 * t62 + t30) * t67 + (-t12 + t14) * qJD(4);];
tauc_reg = t1;
