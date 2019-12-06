% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:28
% EndTime: 2019-12-05 18:16:31
% DurationCPUTime: 0.65s
% Computational Cost: add. (857->119), mult. (1673->171), div. (0->0), fcn. (1058->8), ass. (0->101)
t133 = pkin(7) + pkin(8);
t78 = sin(qJ(5));
t79 = sin(qJ(4));
t81 = cos(qJ(5));
t82 = cos(qJ(4));
t50 = t78 * t82 + t81 * t79;
t73 = qJD(1) + qJD(3);
t45 = t50 * t73;
t72 = qJD(4) + qJD(5);
t131 = qJD(5) - t72;
t129 = pkin(1) * sin(pkin(9));
t103 = qJD(1) * t129;
t66 = cos(pkin(9)) * pkin(1) + pkin(2);
t58 = t66 * qJD(1);
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t39 = t83 * t103 + t80 * t58;
t99 = t133 * t73 + t39;
t19 = t82 * qJD(2) - t99 * t79;
t132 = -t80 * t129 + t83 * t66;
t20 = t79 * qJD(2) + t99 * t82;
t128 = t73 * pkin(3);
t127 = t82 * pkin(4);
t86 = t83 * t129 + t80 * t66;
t47 = pkin(7) + t86;
t126 = -pkin(8) - t47;
t38 = -t80 * t103 + t83 * t58;
t67 = -pkin(3) - t127;
t21 = t67 * t73 - t38;
t108 = t79 * qJD(4);
t102 = pkin(4) * t108;
t110 = qJD(3) * t58;
t95 = qJD(3) * t103;
t36 = t80 * t110 + t83 * t95;
t24 = t73 * t102 + t36;
t28 = t72 * t50;
t115 = t81 * t82;
t119 = t78 * t79;
t49 = -t115 + t119;
t125 = t21 * t28 + t24 * t49;
t27 = t72 * t49;
t124 = -t21 * t27 + t24 * t50;
t22 = t27 * t72;
t123 = t38 * t72;
t122 = t39 * t73;
t41 = t86 * qJD(3);
t121 = t41 * t73;
t104 = t73 * t115;
t105 = t73 * t119;
t43 = -t104 + t105;
t120 = t45 * t43;
t117 = t81 * t20;
t84 = qJD(4) ^ 2;
t113 = t84 * t79;
t107 = t82 * qJD(4);
t30 = -t38 - t128;
t112 = t30 * t107 + t36 * t79;
t111 = t79 ^ 2 - t82 ^ 2;
t106 = pkin(4) * t73 * t79;
t101 = t73 * t107;
t16 = qJD(4) * pkin(4) + t19;
t100 = -pkin(4) * t72 - t16;
t98 = qJD(4) * t133;
t35 = t83 * t110 - t80 * t95;
t97 = -t30 * t73 - t35;
t96 = qJD(4) * t126;
t46 = -pkin(3) - t132;
t94 = -t39 + t102;
t92 = pkin(7) * t84 - t122;
t90 = t47 * t84 + t121;
t89 = qJD(4) * (t38 - t128);
t40 = t132 * qJD(3);
t88 = qJD(4) * (t46 * t73 - t40);
t4 = t19 * qJD(4) + t82 * t35;
t5 = -t20 * qJD(4) - t79 * t35;
t87 = -t21 * t45 - t78 * t4 + t81 * t5;
t17 = qJD(5) * t104 + t81 * t101 - t72 * t105;
t85 = t21 * t43 + (t131 * t20 - t5) * t78;
t71 = t73 ^ 2;
t70 = t82 * pkin(8);
t69 = t84 * t82;
t63 = t82 * pkin(7) + t70;
t62 = t133 * t79;
t53 = 0.2e1 * t79 * t101;
t52 = t82 * t98;
t51 = t79 * t98;
t42 = -0.2e1 * t111 * t73 * qJD(4);
t37 = t46 - t127;
t34 = t41 + t102;
t33 = t82 * t47 + t70;
t32 = t126 * t79;
t25 = t30 * t108;
t23 = t28 * t72;
t18 = t28 * t73;
t12 = -t79 * t40 + t82 * t96;
t11 = t82 * t40 + t79 * t96;
t10 = -t43 ^ 2 + t45 ^ 2;
t6 = t43 * t72 + t17;
t2 = t17 * t50 - t45 * t27;
t1 = -t17 * t49 - t50 * t18 + t27 * t43 - t45 * t28;
t3 = [0, 0, 0, 0, 0, -t36 - t121, -t40 * t73 - t35, t53, t42, t69, -t113, 0, t25 + t79 * t88 + (-t36 - t90) * t82, t90 * t79 + t82 * t88 + t112, t2, t1, -t22, -t23, 0, t34 * t43 + t37 * t18 + (-t78 * t11 + t81 * t12 + (-t32 * t78 - t33 * t81) * qJD(5)) * t72 + t125, t34 * t45 + t37 * t17 - (t81 * t11 + t78 * t12 + (t32 * t81 - t33 * t78) * qJD(5)) * t72 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t69, 0, 0, 0, 0, 0, -t23, t22; 0, 0, 0, 0, 0, -t36 + t122, t38 * t73 - t35, t53, t42, t69, -t113, 0, t25 + t79 * t89 + (-t36 - t92) * t82, t92 * t79 + t82 * t89 + t112, t2, t1, -t22, -t23, 0, t67 * t18 + (t78 * t51 - t81 * t52 + (t62 * t78 - t63 * t81) * qJD(5)) * t72 + t94 * t43 + t50 * t123 + t125, t67 * t17 - (-t81 * t51 - t78 * t52 + (-t62 * t81 - t63 * t78) * qJD(5)) * t72 + t94 * t45 - t49 * t123 + t124; 0, 0, 0, 0, 0, 0, 0, -t79 * t71 * t82, t111 * t71, 0, 0, 0, t97 * t79, t97 * t82, t120, t10, t6, 0, 0, -t43 * t106 - (-t78 * t19 - t117) * t72 + (t100 * t78 - t117) * qJD(5) + t87, -t45 * t106 + (t100 * qJD(5) + t19 * t72 - t4) * t81 + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t10, t6, 0, 0, t87 + t131 * (-t78 * t16 - t117), (-t131 * t16 - t4) * t81 + t85;];
tauc_reg = t3;
