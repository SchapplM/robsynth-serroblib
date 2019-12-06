% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:45
% EndTime: 2019-12-05 16:00:50
% DurationCPUTime: 0.70s
% Computational Cost: add. (1160->159), mult. (2423->222), div. (0->0), fcn. (1530->6), ass. (0->110)
t65 = cos(qJ(2));
t63 = cos(qJ(5));
t64 = cos(qJ(4));
t113 = t63 * t64;
t60 = sin(qJ(5));
t61 = sin(qJ(4));
t76 = t60 * t61 - t113;
t30 = t76 * t65;
t101 = qJD(4) * t64;
t66 = -pkin(2) - pkin(6);
t124 = t66 * qJD(2);
t97 = t65 * qJD(1);
t81 = qJD(3) - t97;
t40 = t81 + t124;
t62 = sin(qJ(2));
t98 = t62 * qJD(1);
t19 = t40 * t101 + (-pkin(7) * t101 + t61 * t98) * qJD(2);
t104 = qJD(2) * t64;
t27 = -pkin(7) * t104 + t64 * t40;
t25 = qJD(4) * pkin(4) + t27;
t126 = (qJD(5) * t25 + t19) * t63;
t58 = t61 ^ 2;
t59 = t64 ^ 2;
t108 = t58 + t59;
t125 = t108 * t62;
t41 = t60 * t64 + t63 * t61;
t34 = t41 * qJD(2);
t57 = qJD(4) + qJD(5);
t123 = pkin(7) - t66;
t43 = t123 * t61;
t44 = t123 * t64;
t23 = -t63 * t43 - t60 * t44;
t102 = qJD(4) * t61;
t38 = t123 * t102;
t39 = qJD(4) * t44;
t72 = t76 * t62;
t122 = -qJD(1) * t72 + t23 * qJD(5) - t63 * t38 - t60 * t39;
t22 = t60 * t43 - t63 * t44;
t73 = t41 * t62;
t121 = qJD(1) * t73 - t22 * qJD(5) - t60 * t38 + t63 * t39;
t100 = qJD(5) * t60;
t79 = t57 * t113;
t21 = -t61 * t100 - t60 * t102 + t79;
t120 = t21 * t57;
t91 = t63 * t104;
t105 = qJD(2) * t61;
t92 = t60 * t105;
t36 = t91 - t92;
t119 = t36 * t34;
t55 = t61 * pkin(4) + qJ(3);
t37 = t55 * qJD(2) + t98;
t118 = t37 * t36;
t45 = (qJD(3) + t97) * qJD(2);
t117 = t45 * t62;
t96 = qJD(2) * qJ(3);
t47 = t96 + t98;
t116 = t47 * t65;
t26 = -pkin(7) * t105 + t61 * t40;
t115 = t60 * t26;
t114 = t63 * t26;
t68 = qJD(2) ^ 2;
t112 = t68 * t62;
t111 = t68 * t65;
t95 = qJD(2) * qJD(4);
t88 = t61 * t95;
t110 = -qJD(5) * t92 - t60 * t88;
t109 = t58 - t59;
t67 = qJD(4) ^ 2;
t107 = t67 + t68;
t106 = qJD(2) * pkin(2);
t103 = qJD(2) * t65;
t99 = t37 * qJD(2);
t94 = t64 * t68 * t61;
t93 = pkin(4) * t104;
t90 = -pkin(4) * t57 - t25;
t89 = t107 * t65;
t87 = t64 * t95;
t49 = t98 * t104;
t18 = t49 + (pkin(7) * qJD(2) - t40) * t102;
t86 = t63 * t18 - t60 * t19;
t85 = -t26 * t100 + t60 * t18;
t84 = -t47 + t98;
t51 = pkin(4) * t101 + qJD(3);
t83 = t51 - t97;
t80 = t61 * t87;
t20 = t57 * t41;
t15 = t20 * qJD(2);
t78 = t15 * t76 - t36 * t20;
t16 = qJD(2) * t79 + t110;
t77 = t41 * t16 + t21 * t34;
t10 = t60 * t25 + t114;
t75 = t45 * qJ(3) + t47 * qJD(3);
t74 = t37 * t34 - t85;
t71 = -t84 + t96;
t1 = t85 + t126;
t2 = -t10 * qJD(5) + t86;
t9 = t63 * t25 - t115;
t70 = -t1 * t41 - t10 * t21 + t2 * t76 + t9 * t20;
t69 = t81 * qJD(2) - t66 * t67 + t45;
t46 = t81 - t106;
t32 = (t51 + t97) * qJD(2);
t31 = t41 * t65;
t17 = t20 * t57;
t13 = -t34 ^ 2 + t36 ^ 2;
t12 = t63 * t27 - t115;
t11 = -t60 * t27 - t114;
t8 = -t110 + (t36 - t91) * t57;
t6 = -qJD(2) * t72 + t20 * t65;
t5 = qJD(2) * t73 + t57 * t30;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t111, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t111, t117 + (t116 + (t46 - t97) * t62) * qJD(2), 0, 0, 0, 0, 0, 0, t61 * t89 + 0.2e1 * t62 * t87, -0.2e1 * t62 * t88 + t64 * t89, -t108 * t112, t117 + (t116 + (t40 - t97) * t125) * qJD(2), 0, 0, 0, 0, 0, 0, t34 * t103 + t62 * t16 + t6 * t57, t36 * t103 - t62 * t15 - t5 * t57, t30 * t15 + t31 * t16 - t5 * t34 - t6 * t36, -t1 * t31 + t10 * t5 + t2 * t30 + t32 * t62 + t9 * t6 + t65 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), (-t116 + (-t46 - t106) * t62) * qJD(1) + t75, -0.2e1 * t80, 0.2e1 * t109 * t95, -t67 * t61, 0.2e1 * t80, -t67 * t64, 0, t71 * t101 + t69 * t61, -t71 * t102 + t69 * t64, 0, (-t116 + (-t40 + t124) * t125) * qJD(1) + t75, t78, t15 * t41 + t16 * t76 + t20 * t34 - t36 * t21, -t17, t77, -t120, 0, -t122 * t57 + t55 * t16 + t37 * t21 + t32 * t41 + t83 * t34, t121 * t57 - t55 * t15 - t37 * t20 - t32 * t76 + t83 * t36, t121 * t34 + t122 * t36 + t22 * t15 - t23 * t16 + t70, t1 * t23 - t121 * t10 - t122 * t9 + t2 * t22 + t32 * t55 + t83 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t84 * qJD(2), 0, 0, 0, 0, 0, 0, -t107 * t61, -t107 * t64, 0, (t108 * t98 - t47) * qJD(2), 0, 0, 0, 0, 0, 0, -qJD(2) * t34 - t17, -qJD(2) * t36 - t120, -t77 - t78, -t70 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t109 * t68, 0, -t94, 0, 0, -t47 * t104 + t49, -t84 * t105, 0, 0, t119, t13, 0, -t119, t8, 0, -t34 * t93 - t11 * t57 - t118 + (t90 * t60 - t114) * qJD(5) + t86, -t36 * t93 + t12 * t57 + (t90 * qJD(5) - t19) * t63 + t74, (t10 + t11) * t36 + (t12 - t9) * t34 + (t15 * t63 - t16 * t60 + (-t34 * t63 + t36 * t60) * qJD(5)) * pkin(4), -t10 * t12 - t9 * t11 + (-t64 * t99 + t1 * t60 + t2 * t63 + (t10 * t63 - t60 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t13, 0, -t119, t8, 0, t10 * t57 - t118 + t2, t9 * t57 - t126 + t74, 0, 0;];
tauc_reg = t3;
