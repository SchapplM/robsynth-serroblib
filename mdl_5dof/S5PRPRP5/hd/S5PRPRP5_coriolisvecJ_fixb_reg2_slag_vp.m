% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:41
% EndTime: 2019-12-05 15:38:44
% DurationCPUTime: 0.70s
% Computational Cost: add. (1034->153), mult. (2740->193), div. (0->0), fcn. (1978->6), ass. (0->98)
t69 = cos(qJ(2));
t65 = sin(pkin(8));
t67 = sin(qJ(4));
t115 = t67 * t65;
t118 = cos(qJ(4));
t66 = cos(pkin(8));
t97 = t118 * t66;
t80 = t97 - t115;
t77 = t80 * t69;
t112 = pkin(6) + qJ(3);
t53 = t112 * t65;
t54 = t112 * t66;
t81 = -t118 * t53 - t67 * t54;
t110 = -qJD(1) * t77 + qJD(3) * t80 + qJD(4) * t81;
t50 = t118 * t65 + t67 * t66;
t101 = t69 * qJD(1);
t51 = (qJD(3) + t101) * qJD(2);
t129 = t50 * t51;
t76 = qJD(2) * t50;
t86 = qJD(2) * t97;
t98 = qJD(2) * t115;
t41 = -t86 + t98;
t61 = -t66 * pkin(3) - pkin(2);
t105 = qJD(2) * t61;
t88 = qJD(3) - t101;
t47 = t88 + t105;
t9 = t41 * pkin(4) - qJ(5) * t76 + t47;
t130 = t9 * t76 + t129;
t93 = qJD(4) * t118;
t127 = -qJD(4) * t115 + t66 * t93;
t107 = t65 ^ 2 + t66 ^ 2;
t68 = sin(qJ(2));
t102 = t68 * qJD(1);
t57 = qJD(2) * qJ(3) + t102;
t126 = t107 * t57;
t122 = t76 ^ 2;
t37 = t41 ^ 2;
t125 = -t37 - t122;
t124 = -t37 + t122;
t46 = t50 * qJD(4);
t123 = t110 * qJD(4);
t92 = pkin(6) * qJD(2) + t57;
t30 = t92 * t65;
t31 = t92 * t66;
t16 = t118 * t31 - t67 * t30;
t3 = t16 * qJD(4) + t129;
t121 = t3 * t81;
t34 = t50 * t68;
t120 = t3 * t34;
t117 = t76 * t41;
t116 = t67 * t31;
t70 = qJD(2) ^ 2;
t114 = t70 * t68;
t113 = t70 * t69;
t11 = qJD(4) * qJ(5) + t16;
t111 = t11 - t16;
t25 = t118 * t54 - t67 * t53;
t109 = qJD(4) * t25 + t88 * t50;
t108 = -t30 * t93 + t51 * t97;
t106 = qJD(2) * pkin(2);
t104 = qJD(2) * t68;
t103 = qJD(4) * t46;
t15 = -t118 * t30 - t116;
t100 = qJD(5) - t15;
t96 = t107 * t69;
t95 = t107 * t51;
t94 = t15 + t116;
t91 = -t46 * pkin(4) + qJ(5) * t127 + t50 * qJD(5) + t102;
t90 = t109 * qJD(4);
t89 = t107 * qJD(3);
t33 = qJD(2) * t46;
t85 = -t33 * t80 + t41 * t46;
t84 = t51 * t115 - t108;
t56 = qJD(4) * t86;
t32 = qJD(4) * t98 - t56;
t62 = qJD(2) * t102;
t83 = t33 * pkin(4) + t32 * qJ(5) + t62;
t18 = t127 * t68 + t69 * t76;
t82 = -t18 * qJD(4) + t41 * t104 - t69 * t33;
t17 = qJD(2) * t77 - t68 * t46;
t35 = t80 * t68;
t79 = -t17 * t41 + t18 * t76 - t34 * t32 - t35 * t33;
t75 = t17 * qJD(4) - t104 * t76 - t69 * t32;
t74 = -t127 * t41 - t32 * t80 - t50 * t33 - t46 * t76;
t73 = 0.2e1 * t76 * qJD(4);
t72 = t109 * t76 - t110 * t41 - t25 * t33 + t3 * t50 + t81 * t32;
t55 = t88 - t106;
t36 = t127 * qJD(4);
t23 = -pkin(4) * t80 - t50 * qJ(5) + t61;
t22 = pkin(4) * t76 + t41 * qJ(5);
t21 = t56 + (t41 - t98) * qJD(4);
t20 = -t56 + (t41 + t98) * qJD(4);
t10 = -qJD(4) * pkin(4) + t100;
t5 = -qJD(5) * t76 + t83;
t4 = t127 * t76 - t32 * t50;
t2 = (-qJD(4) * t31 - t51 * t65) * t67 + t108;
t1 = (qJD(5) - t116) * qJD(4) - t84;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -t113, 0, 0, 0, 0, 0, 0, 0, 0, -t66 * t114, t65 * t114, t107 * t113, t68 * t95 + (t55 * t68 + (-t102 + t126) * t69) * qJD(2), 0, 0, 0, 0, 0, 0, t82, -t75, t79, -t15 * t18 + t16 * t17 + t2 * t35 + t120 + (t47 - t101) * t104, 0, 0, 0, 0, 0, 0, t82, t79, t75, t1 * t35 + t10 * t18 + t9 * t104 + t11 * t17 - t5 * t69 + t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 + (-qJD(1) * t96 + t89) * qJD(2), t57 * t89 + qJ(3) * t95 + ((-t55 - t106) * t68 - t57 * t96) * qJD(1), t4, t74, t36, t85, -t103, 0, t61 * t33 + t47 * t46 - t90 + (-qJD(2) * t80 - t41) * t102, t127 * t47 - t61 * t32 - t123, -t127 * t15 - t16 * t46 + t2 * t80 + t72, t2 * t25 - t121 + t110 * t16 - t109 * t15 + (-t47 + t105) * t102, t4, t36, -t74, 0, t103, t85, t23 * t33 - t91 * t41 + t9 * t46 - t5 * t80 - t90, t1 * t80 + t10 * t127 - t11 * t46 + t72, -t127 * t9 + t23 * t32 - t5 * t50 + t76 * t91 + t123, t1 * t25 + t109 * t10 + t110 * t11 + t5 * t23 - t91 * t9 - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107 * t70, -qJD(2) * t126 + t62, 0, 0, 0, 0, 0, 0, t73, -t20, t125, t15 * t76 + t16 * t41 + t62, 0, 0, 0, 0, 0, 0, t73, t125, t20, t11 * t41 + (-qJD(5) - t10) * t76 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, t124, t21, -t117, 0, 0, -t47 * t76 - t129, t94 * qJD(4) + t47 * t41 + t84, 0, 0, t117, t21, -t124, 0, 0, -t117, -t22 * t41 - t130, pkin(4) * t32 - t33 * qJ(5) + t111 * t76 + (t10 - t100) * t41, t22 * t76 - t9 * t41 + (0.2e1 * qJD(5) - t94) * qJD(4) - t84, -t3 * pkin(4) + t1 * qJ(5) - t10 * t16 + t100 * t11 - t9 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, t21, -qJD(4) ^ 2 - t122, -t111 * qJD(4) + t130;];
tauc_reg = t6;
