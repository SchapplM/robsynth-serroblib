% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:54
% EndTime: 2019-12-05 15:57:57
% DurationCPUTime: 0.90s
% Computational Cost: add. (1042->175), mult. (2883->269), div. (0->0), fcn. (2327->10), ass. (0->106)
t69 = cos(pkin(10));
t106 = qJD(2) * t69;
t75 = cos(qJ(4));
t60 = t75 * t106;
t67 = sin(pkin(10));
t72 = sin(qJ(4));
t114 = t72 * t67;
t96 = qJD(2) * t114;
t43 = -t60 + t96;
t40 = qJD(5) + t43;
t68 = sin(pkin(5));
t108 = qJD(1) * t68;
t76 = cos(qJ(2));
t97 = t76 * t108;
t49 = (qJD(3) + t97) * qJD(2);
t50 = -t75 * t69 + t114;
t130 = t50 * t49;
t51 = t75 * t67 + t72 * t69;
t129 = t51 * t49;
t89 = qJD(3) - t97;
t128 = -qJD(5) + t40;
t73 = sin(qJ(2));
t98 = t73 * t108;
t53 = qJD(2) * qJ(3) + t98;
t70 = cos(pkin(5));
t107 = qJD(1) * t70;
t59 = t69 * t107;
t24 = t59 + (-pkin(7) * qJD(2) - t53) * t67;
t34 = t67 * t107 + t69 * t53;
t25 = pkin(7) * t106 + t34;
t10 = t72 * t24 + t75 * t25;
t4 = t10 * qJD(4) + t129;
t45 = qJD(2) * t51;
t127 = (t45 * pkin(4) + t40 * pkin(8)) * t40 + t4;
t62 = -t69 * pkin(3) - pkin(2);
t39 = t62 * qJD(2) + t89;
t11 = t43 * pkin(4) - t45 * pkin(8) + t39;
t118 = t68 * t76;
t80 = t50 * t118;
t113 = pkin(7) + qJ(3);
t54 = t113 * t67;
t55 = t113 * t69;
t86 = -t75 * t54 - t72 * t55;
t112 = -qJD(1) * t80 + t50 * qJD(3) - t86 * qJD(4);
t23 = t50 * pkin(4) - t51 * pkin(8) + t62;
t27 = -t72 * t54 + t75 * t55;
t9 = t75 * t24 - t72 * t25;
t3 = t9 * qJD(4) - t130;
t47 = t51 * qJD(4);
t38 = qJD(2) * t47;
t46 = t50 * qJD(4);
t5 = -qJD(4) * pkin(4) - t9;
t126 = -(qJD(5) * t11 + t3) * t50 + t4 * t51 - t5 * t46 + (-qJD(5) * t23 + t112) * t40 - t27 * t38;
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t32 = t71 * qJD(4) + t74 * t45;
t57 = qJD(4) * t60;
t37 = -qJD(4) * t96 + t57;
t15 = t32 * qJD(5) + t71 * t37;
t102 = t74 * qJD(4);
t103 = qJD(5) * t71;
t14 = qJD(5) * t102 - t45 * t103 + t74 * t37;
t125 = t14 * t71;
t124 = t23 * t38;
t30 = t71 * t45 - t102;
t123 = t30 * t40;
t122 = t32 * t40;
t121 = t32 * t45;
t120 = t45 * t30;
t119 = t68 * t73;
t77 = qJD(2) ^ 2;
t117 = t68 * t77;
t115 = t71 * t38;
t36 = t74 * t38;
t81 = t51 * t118;
t111 = -qJD(1) * t81 + t51 * qJD(3) + t27 * qJD(4);
t110 = t67 ^ 2 + t69 ^ 2;
t109 = qJD(2) * pkin(2);
t105 = qJD(2) * t73;
t104 = qJD(5) * t51;
t101 = t73 * t117;
t100 = t76 * t117;
t99 = t68 * t105;
t95 = t110 * t49;
t92 = t40 * t74;
t90 = t47 * pkin(4) + t46 * pkin(8) - t98;
t6 = qJD(4) * pkin(8) + t10;
t1 = t74 * t11 - t71 * t6;
t2 = t71 * t11 + t74 * t6;
t88 = (-t67 * t53 + t59) * t67 - t34 * t69;
t41 = -t67 * t119 + t70 * t69;
t42 = t69 * t119 + t70 * t67;
t87 = t75 * t41 - t72 * t42;
t19 = t72 * t41 + t75 * t42;
t85 = t36 + (-t43 * t71 - t103) * t40;
t84 = -t74 * t118 - t71 * t19;
t83 = t71 * t118 - t74 * t19;
t82 = -t51 * t103 - t74 * t46;
t79 = -pkin(8) * t38 + (t5 + t9) * t40;
t56 = qJD(2) * t98;
t52 = t89 - t109;
t13 = t38 * pkin(4) - t37 * pkin(8) + t56;
t12 = t74 * t13;
t8 = qJD(2) * t81 + t19 * qJD(4);
t7 = -qJD(2) * t80 + t87 * qJD(4);
t16 = [0, 0, -t101, -t100, -t69 * t101, t67 * t101, t110 * t100, (-t41 * t67 + t42 * t69) * t49 + (t52 * t73 + (-t88 - t98) * t76) * t68 * qJD(2), 0, 0, 0, 0, 0, -t8 * qJD(4) + (t43 * t105 - t38 * t76) * t68, -t7 * qJD(4) + (t45 * t105 - t37 * t76) * t68, 0, 0, 0, 0, 0, (t83 * qJD(5) - t71 * t7 + t74 * t99) * t40 + t84 * t38 + t8 * t30 - t87 * t15, -(t84 * qJD(5) + t74 * t7 + t71 * t99) * t40 + t83 * t38 + t8 * t32 - t87 * t14; 0, 0, 0, 0, 0, 0, qJD(2) * t110 * t89 + t95, -t88 * qJD(3) + qJ(3) * t95 + (t88 * t76 + (-t52 - t109) * t73) * t108, t37 * t51 - t45 * t46, -t37 * t50 - t51 * t38 + t46 * t43 - t45 * t47, -t46 * qJD(4), -t47 * qJD(4), 0, t62 * t38 + t39 * t47 - t111 * qJD(4) + (qJD(2) * t50 - t43) * t98, t112 * qJD(4) + t62 * t37 - t39 * t46, t14 * t74 * t51 + t82 * t32, -(-t30 * t74 - t32 * t71) * t46 + (-t125 - t15 * t74 + (t30 * t71 - t32 * t74) * qJD(5)) * t51, t14 * t50 + t32 * t47 + t51 * t36 + t82 * t40, -t51 * t115 - t15 * t50 - t30 * t47 + (-t74 * t104 + t71 * t46) * t40, t38 * t50 + t40 * t47, t1 * t47 + t12 * t50 - t86 * t15 + t111 * t30 + (t124 + t90 * t40 + (-t27 * t40 + t5 * t51 - t6 * t50) * qJD(5)) * t74 + t126 * t71, -t86 * t14 - t2 * t47 + t111 * t32 + (-t124 - (-qJD(5) * t6 + t13) * t50 - t5 * t104 + (qJD(5) * t27 - t90) * t40) * t71 + t126 * t74; 0, 0, 0, 0, 0, 0, -t110 * t77, t88 * qJD(2) + t56, 0, 0, 0, 0, 0, 0.2e1 * t45 * qJD(4), t57 + (-t43 - t96) * qJD(4), 0, 0, 0, 0, 0, t85 - t120, -t40 ^ 2 * t74 - t115 - t121; 0, 0, 0, 0, 0, 0, 0, 0, t45 * t43, -t43 ^ 2 + t45 ^ 2, t57 + (t43 - t96) * qJD(4), 0, 0, -t39 * t45 - t129, t39 * t43 + t130, t32 * t92 + t125, (t14 - t123) * t74 + (-t15 - t122) * t71, t40 * t92 + t115 - t121, t85 + t120, -t40 * t45, -pkin(4) * t15 - t1 * t45 - t10 * t30 - t127 * t74 + t79 * t71, -pkin(4) * t14 - t10 * t32 + t127 * t71 + t2 * t45 + t79 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t30, -t30 ^ 2 + t32 ^ 2, t14 + t123, t122 - t15, t38, t128 * t2 - t71 * t3 - t5 * t32 + t12, t1 * t128 - t71 * t13 - t74 * t3 + t5 * t30;];
tauc_reg = t16;
