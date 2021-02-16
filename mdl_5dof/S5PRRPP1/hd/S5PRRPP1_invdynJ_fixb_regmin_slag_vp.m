% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:22:40
% EndTime: 2021-01-15 15:22:44
% DurationCPUTime: 0.84s
% Computational Cost: add. (1276->214), mult. (2679->265), div. (0->0), fcn. (1753->8), ass. (0->106)
t78 = pkin(7) + qJ(2);
t72 = sin(t78);
t74 = cos(t78);
t141 = g(1) * t72 - g(2) * t74;
t109 = g(1) * t74 + g(2) * t72;
t83 = sin(pkin(8));
t84 = cos(pkin(8));
t86 = sin(qJ(3));
t87 = cos(qJ(3));
t52 = t83 * t87 + t84 * t86;
t48 = t52 * qJD(2);
t44 = t48 ^ 2;
t130 = t84 * t87;
t116 = qJD(2) * t130;
t127 = qJD(2) * t86;
t45 = t83 * t127 - t116;
t142 = -t45 ^ 2 - t44;
t129 = qJ(4) + pkin(6);
t114 = t129 * t86;
t61 = t129 * t87;
t33 = -t83 * t114 + t84 * t61;
t79 = qJ(3) + pkin(8);
t73 = sin(t79);
t140 = -t33 * qJDD(3) - t141 * t73;
t111 = t129 * qJD(3);
t106 = qJD(2) * t111;
t138 = qJD(3) * qJD(1) + qJD(2) * qJD(4) + t129 * qJDD(2);
t76 = t87 * qJDD(1);
t17 = qJDD(3) * pkin(3) - t87 * t106 - t138 * t86 + t76;
t21 = (qJDD(1) - t106) * t86 + t138 * t87;
t5 = t83 * t17 + t84 * t21;
t75 = cos(t79);
t139 = g(3) * t73 + t109 * t75 - t5;
t137 = pkin(3) * t86;
t133 = g(3) * t87;
t132 = t87 * pkin(3);
t42 = t86 * qJD(1) + qJD(2) * t61;
t131 = t83 * t42;
t35 = t84 * t42;
t4 = t84 * t17 - t83 * t21;
t40 = t87 * qJD(1) - qJD(2) * t114;
t38 = qJD(3) * pkin(3) + t40;
t20 = t83 * t38 + t35;
t81 = t86 ^ 2;
t128 = -t87 ^ 2 + t81;
t126 = qJDD(3) * pkin(4);
t23 = t83 * t40 + t35;
t125 = t23 * qJD(3);
t124 = t86 * qJD(3);
t25 = t84 * t40 - t131;
t123 = qJD(5) - t25;
t122 = qJDD(1) - g(3);
t120 = t86 * qJDD(2);
t119 = t87 * qJDD(2);
t118 = qJD(2) * qJD(3);
t117 = pkin(3) * t124;
t71 = pkin(2) + t132;
t113 = t86 * t118;
t112 = t87 * t118;
t107 = -t84 * t119 + t83 * t120;
t105 = t75 * pkin(4) + t73 * qJ(5);
t19 = t84 * t38 - t131;
t47 = t52 * qJD(3);
t51 = t83 * t86 - t130;
t104 = -t47 * qJD(3) - t51 * qJDD(3);
t103 = -g(3) * t75 + t109 * t73 + t4;
t101 = -0.2e1 * pkin(2) * t118 - pkin(6) * qJDD(3);
t58 = -t71 * qJD(2) + qJD(4);
t28 = qJD(2) * t47 + t107;
t96 = t52 * qJDD(2) - t83 * t113;
t29 = t84 * t112 + t96;
t50 = qJD(3) * t130 - t83 * t124;
t100 = -t52 * t28 + t51 * t29 - t50 * t45 + t47 * t48;
t32 = t84 * t114 + t83 * t61;
t99 = -t32 * qJDD(3) + t141 * t75;
t97 = -t86 * qJD(4) - t87 * t111;
t39 = pkin(3) * t113 - t71 * qJDD(2) + qJDD(4);
t88 = qJD(3) ^ 2;
t95 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t88 + t141;
t89 = qJD(2) ^ 2;
t94 = t89 * pkin(2) - qJDD(2) * pkin(6) + t109;
t18 = t45 * pkin(4) - t48 * qJ(5) + t58;
t93 = -t18 * t48 - qJDD(5) + t103;
t41 = t87 * qJD(4) - t86 * t111;
t24 = t83 * t41 - t84 * t97;
t26 = t84 * t41 + t83 * t97;
t92 = t24 * t48 - t26 * t45 - t33 * t28 + t32 * t29 - t109;
t91 = t28 * pkin(4) - t29 * qJ(5) + t39;
t90 = 0.2e1 * t48 * qJD(3) + t107;
t80 = qJDD(3) * qJ(5);
t70 = -t84 * pkin(3) - pkin(4);
t67 = t83 * pkin(3) + qJ(5);
t60 = qJDD(3) * t87 - t88 * t86;
t59 = qJDD(3) * t86 + t88 * t87;
t57 = t74 * t71;
t30 = t50 * qJD(3) + t52 * qJDD(3);
t27 = t51 * pkin(4) - t52 * qJ(5) - t71;
t22 = pkin(3) * t127 + t48 * pkin(4) + t45 * qJ(5);
t14 = qJD(3) * qJ(5) + t20;
t12 = -qJD(3) * pkin(4) + qJD(5) - t19;
t9 = (-t45 + t116) * qJD(3) + t96;
t8 = t47 * pkin(4) - t50 * qJ(5) - t52 * qJD(5) + t117;
t3 = qJDD(5) - t126 - t4;
t2 = qJD(3) * qJD(5) + t5 + t80;
t1 = -t48 * qJD(5) + t91;
t6 = [t122, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t59, t104, -t30, t100, -t19 * t47 + t20 * t50 - t4 * t51 + t5 * t52 - g(3), t104, t100, t30, t12 * t47 + t14 * t50 + t2 * t52 + t3 * t51 - g(3); 0, qJDD(2), t141, t109, t81 * qJDD(2) + 0.2e1 * t86 * t112, -0.2e1 * t128 * t118 + 0.2e1 * t86 * t119, t59, t60, 0, t101 * t86 + t95 * t87, t101 * t87 - t95 * t86, -t71 * t28 + t39 * t51 + t58 * t47 + (t45 * t137 - t24) * qJD(3) + t99, -t71 * t29 + t39 * t52 + t58 * t50 + (t48 * t137 - t26) * qJD(3) + t140, -t19 * t50 - t20 * t47 - t4 * t52 - t5 * t51 + t92, t5 * t33 + t20 * t26 - t4 * t32 - t19 * t24 - t39 * t71 + t58 * t117 - g(1) * (t129 * t74 - t72 * t71) - g(2) * (t129 * t72 + t57), -t24 * qJD(3) + t1 * t51 + t18 * t47 + t27 * t28 + t8 * t45 + t99, t12 * t50 - t14 * t47 - t2 * t51 + t3 * t52 + t92, t26 * qJD(3) - t1 * t52 - t18 * t50 - t27 * t29 - t8 * t48 - t140, -g(2) * t57 + t1 * t27 + t12 * t24 + t14 * t26 + t18 * t8 + t2 * t33 + t3 * t32 + (-g(1) * t129 - g(2) * t105) * t74 + (-g(1) * (-t105 - t71) - g(2) * t129) * t72; 0, 0, 0, 0, -t86 * t89 * t87, t128 * t89, t120, t119, qJDD(3), t86 * t94 - t133 + t76, -t122 * t86 + t94 * t87, t125 - t58 * t48 + (qJDD(3) * t84 - t45 * t127) * pkin(3) + t103, t25 * qJD(3) + t58 * t45 + (-qJDD(3) * t83 - t48 * t127) * pkin(3) + t139, (t20 - t23) * t48 + (-t19 + t25) * t45 + (-t28 * t83 - t29 * t84) * pkin(3), t19 * t23 - t20 * t25 + (-t133 + t4 * t84 + t5 * t83 + (-qJD(2) * t58 + t109) * t86) * pkin(3), t125 - t22 * t45 + (pkin(4) - t70) * qJDD(3) + t93, -t67 * t28 + t70 * t29 + (t14 - t23) * t48 + (t12 - t123) * t45, t67 * qJDD(3) - t18 * t45 + t22 * t48 + t80 + (0.2e1 * qJD(5) - t25) * qJD(3) - t139, t2 * t67 + t3 * t70 - t18 * t22 - t12 * t23 - g(3) * (t105 + t132) + t123 * t14 + t109 * (pkin(4) * t73 - qJ(5) * t75 + t137); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t9, t142, t19 * t48 + t20 * t45 - t141 + t39, t90, t142, -t9, t14 * t45 + (-qJD(5) - t12) * t48 + t91 - t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t45 - qJDD(3), (t45 + t116) * qJD(3) + t96, -t44 - t88, -t14 * qJD(3) - t126 - t93;];
tau_reg = t6;
