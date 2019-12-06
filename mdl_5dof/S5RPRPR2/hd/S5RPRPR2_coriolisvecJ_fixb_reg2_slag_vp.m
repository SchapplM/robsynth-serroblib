% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:38
% EndTime: 2019-12-05 17:49:41
% DurationCPUTime: 0.67s
% Computational Cost: add. (1513->126), mult. (3016->174), div. (0->0), fcn. (1930->8), ass. (0->93)
t88 = sin(pkin(9));
t90 = cos(pkin(9));
t111 = t88 ^ 2 + t90 ^ 2;
t128 = sin(pkin(8)) * pkin(1);
t107 = qJD(3) * t128;
t101 = qJD(1) * t107;
t95 = cos(qJ(3));
t110 = qJD(3) * t95;
t79 = cos(pkin(8)) * pkin(1) + pkin(2);
t72 = t79 * qJD(1);
t93 = sin(qJ(3));
t113 = t93 * t101 - t72 * t110;
t87 = qJD(1) + qJD(3);
t34 = t87 * qJD(4) - t113;
t134 = t111 * t34;
t108 = qJD(1) * t128;
t46 = -t93 * t108 + t95 * t72;
t99 = qJD(4) - t46;
t94 = cos(qJ(5));
t120 = t94 * t90;
t92 = sin(qJ(5));
t122 = t92 * t88;
t63 = -t120 + t122;
t133 = t63 * t34;
t64 = t94 * t88 + t92 * t90;
t52 = t64 * t87;
t112 = t95 * t128 + t93 * t79;
t132 = t111 * t87;
t121 = t93 * t72;
t42 = qJD(3) * t121 + t95 * t101;
t55 = t112 * qJD(3);
t131 = -t55 * t87 - t42;
t47 = t95 * t108 + t121;
t130 = t47 * t87 - t42;
t129 = t52 ^ 2;
t127 = t90 * pkin(4);
t83 = t90 * pkin(7);
t109 = t87 * t122;
t50 = -t87 * t120 + t109;
t125 = t52 * t50;
t123 = t87 * t88;
t67 = (-pkin(7) - qJ(4)) * t88;
t68 = t90 * qJ(4) + t83;
t39 = t92 * t67 + t94 * t68;
t119 = t39 * qJD(5) + t99 * t64;
t38 = t94 * t67 - t92 * t68;
t118 = -t38 * qJD(5) + t99 * t63;
t80 = -pkin(3) - t127;
t30 = t80 * t87 + t99;
t60 = t64 * qJD(5);
t117 = t30 * t60 + t42 * t63;
t105 = qJD(5) * t120;
t106 = qJD(5) * t122;
t59 = -t105 + t106;
t116 = -t30 * t59 + t42 * t64;
t44 = t87 * t60;
t115 = -t64 * t44 + t59 * t50;
t37 = t87 * qJ(4) + t47;
t27 = t88 * qJD(2) + t90 * t37;
t82 = t90 * qJD(2);
t22 = t82 + (-pkin(7) * t87 - t37) * t88;
t23 = t87 * t83 + t27;
t10 = t94 * t22 - t92 * t23;
t11 = t92 * t22 + t94 * t23;
t3 = t10 * qJD(5) - t133;
t96 = t64 * t34;
t4 = -t11 * qJD(5) - t96;
t104 = t10 * t59 - t11 * t60 - t3 * t63 - t4 * t64;
t102 = -t93 * t128 + t95 * t79;
t100 = -pkin(3) - t102;
t98 = (-t88 * t37 + t82) * t88 - t27 * t90;
t56 = qJ(4) + t112;
t40 = (-pkin(7) - t56) * t88;
t41 = t90 * t56 + t83;
t14 = t94 * t40 - t92 * t41;
t15 = t92 * t40 + t94 * t41;
t65 = t87 * t105;
t43 = t87 * t106 - t65;
t97 = -t63 * t43 + t52 * t60;
t54 = -t93 * t107 + t79 * t110;
t58 = t60 * qJD(5);
t57 = t59 * qJD(5);
t49 = t50 ^ 2;
t48 = qJD(4) + t54;
t45 = t100 - t127;
t36 = t42 * t88;
t35 = -t87 * pkin(3) + t99;
t13 = t44 * t63 + t50 * t60;
t12 = -t43 * t64 - t52 * t59;
t7 = -t15 * qJD(5) - t64 * t48;
t6 = t14 * qJD(5) - t63 * t48;
t5 = -t97 + t115;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, -t54 * t87 + t113, 0, -t42 * t102 - t113 * t112 - t46 * t55 + t47 * t54, 0, 0, 0, 0, 0, 0, t131 * t90, t55 * t123 + t36, t48 * t132 + t134, t42 * t100 + t134 * t56 + t35 * t55 - t98 * t48, t12, t5, -t57, t13, -t58, 0, t7 * qJD(5) + t45 * t44 + t55 * t50 + t117, -t6 * qJD(5) - t45 * t43 + t55 * t52 + t116, t14 * t43 - t15 * t44 - t6 * t50 - t7 * t52 + t104, t10 * t7 + t11 * t6 + t4 * t14 + t3 * t15 + t30 * t55 + t42 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t57, t97 + t115, -t10 * t60 - t11 * t59 + t3 * t64 - t4 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t46 * t87 + t113, 0, 0, 0, 0, 0, 0, 0, 0, t130 * t90, -t47 * t123 + t36, t99 * t132 + t134, -t42 * pkin(3) + qJ(4) * t134 - t35 * t47 - t99 * t98, t12, t5, -t57, t13, -t58, 0, -t119 * qJD(5) + t80 * t44 - t47 * t50 + t117, t118 * qJD(5) - t80 * t43 - t47 * t52 + t116, t118 * t50 + t119 * t52 + t38 * t43 - t39 * t44 + t104, -t119 * t10 - t118 * t11 + t3 * t39 - t30 * t47 + t4 * t38 + t42 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111 * t87 ^ 2, t98 * t87 + t42, 0, 0, 0, 0, 0, 0, 0.2e1 * t52 * qJD(5), t65 + (-t50 - t109) * qJD(5), -t49 - t129, t10 * t52 + t11 * t50 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, -t49 + t129, t65 + (t50 - t109) * qJD(5), -t125, 0, 0, -t30 * t52 - t96, t30 * t50 + t133, 0, 0;];
tauc_reg = t1;
