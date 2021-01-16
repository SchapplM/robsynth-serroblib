% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:52
% EndTime: 2021-01-15 11:33:57
% DurationCPUTime: 0.87s
% Computational Cost: add. (1113->170), mult. (2507->252), div. (0->0), fcn. (1710->6), ass. (0->108)
t103 = sin(qJ(5));
t105 = cos(qJ(5));
t123 = qJD(5) * t103;
t101 = sin(pkin(8));
t102 = cos(pkin(8));
t104 = sin(qJ(3));
t106 = cos(qJ(3));
t78 = t101 * t106 + t102 * t104;
t68 = t78 * qJD(1);
t135 = t105 * t68;
t127 = qJD(1) * t104;
t116 = t101 * t127;
t120 = qJD(1) * qJD(3);
t115 = t106 * t120;
t84 = t102 * t115;
t62 = qJD(3) * t116 - t84;
t63 = t78 * t120;
t126 = qJD(1) * t106;
t72 = t102 * t126 - t116;
t1 = -qJD(5) * t135 + t103 * t62 - t105 * t63 - t72 * t123;
t26 = t103 * t72 + t135;
t96 = qJD(3) + qJD(5);
t137 = t26 * t96;
t152 = t1 + t137;
t112 = -t103 * t68 + t105 * t72;
t151 = t112 * t26;
t138 = t112 * t96;
t2 = t112 * qJD(5) - t103 * t63 - t105 * t62;
t150 = -t2 + t138;
t149 = t112 ^ 2 - t26 ^ 2;
t142 = t68 * pkin(7);
t107 = -pkin(1) - pkin(6);
t85 = t107 * qJD(1) + qJD(2);
t66 = -qJ(4) * t127 + t104 * t85;
t136 = t102 * t66;
t67 = -qJ(4) * t126 + t106 * t85;
t55 = qJD(3) * pkin(3) + t67;
t18 = t101 * t55 + t136;
t10 = t18 - t142;
t81 = pkin(3) * t127 + qJD(1) * qJ(2) + qJD(4);
t38 = t68 * pkin(4) + t81;
t148 = t10 * t123 + t38 * t26;
t146 = qJD(5) - t96;
t121 = t106 * qJD(4);
t125 = qJD(3) * t104;
t39 = -t85 * t125 + (qJ(4) * t125 - t121) * qJD(1);
t122 = t104 * qJD(4);
t124 = qJD(3) * t106;
t40 = t85 * t124 + (-qJ(4) * t124 - t122) * qJD(1);
t13 = -t101 * t40 + t102 * t39;
t4 = t63 * pkin(7) + t13;
t14 = t101 * t39 + t102 * t40;
t5 = t62 * pkin(7) + t14;
t145 = -t103 * t5 + t105 * t4 - t38 * t112;
t97 = qJD(1) * qJD(2);
t144 = 0.2e1 * t97;
t77 = -t101 * t104 + t102 * t106;
t30 = t103 * t77 + t105 * t78;
t70 = t101 * t125 - t102 * t124;
t71 = t78 * qJD(3);
t6 = -t30 * qJD(5) + t103 * t70 - t105 * t71;
t143 = t6 * t96;
t31 = -t103 * t78 + t105 * t77;
t7 = t31 * qJD(5) - t103 * t71 - t105 * t70;
t141 = t7 * t96;
t140 = t72 * pkin(7);
t139 = pkin(3) * t101;
t129 = qJ(4) - t107;
t64 = t129 * t125 - t121;
t83 = t129 * t106;
t65 = -qJD(3) * t83 - t122;
t20 = t101 * t64 + t102 * t65;
t51 = t101 * t66;
t24 = t102 * t67 - t51;
t82 = t129 * t104;
t37 = -t101 * t83 - t102 * t82;
t80 = pkin(3) * t115 + t97;
t134 = -t104 ^ 2 + t106 ^ 2;
t90 = t104 * pkin(3) + qJ(2);
t108 = qJD(3) ^ 2;
t133 = t108 * t104;
t132 = t108 * t106;
t109 = qJD(1) ^ 2;
t131 = t109 * qJ(2);
t130 = t109 * t106;
t128 = -t108 - t109;
t86 = pkin(3) * t124 + qJD(2);
t119 = 0.2e1 * qJD(1);
t118 = pkin(3) * t126;
t17 = t102 * t55 - t51;
t19 = -t101 * t65 + t102 * t64;
t23 = -t101 * t67 - t136;
t36 = t101 * t82 - t102 * t83;
t9 = qJD(3) * pkin(4) - t140 + t17;
t113 = -t105 * t10 - t103 * t9;
t111 = t13 * t77 + t14 * t78 - t17 * t71 - t18 * t70;
t91 = t102 * pkin(3) + pkin(4);
t56 = t78 * pkin(4) + t90;
t44 = t72 * pkin(4) + t118;
t41 = -t70 * pkin(4) + t86;
t32 = -t62 * pkin(4) + t80;
t22 = -t78 * pkin(7) + t37;
t21 = -t77 * pkin(7) + t36;
t16 = t24 - t140;
t15 = t23 + t142;
t12 = t70 * pkin(7) + t20;
t11 = t71 * pkin(7) + t19;
t3 = [0, 0, 0, 0, t144, qJ(2) * t144, -0.2e1 * t104 * t115, -0.2e1 * t134 * t120, -t133, -t132, 0, -t107 * t133 + (qJ(2) * t124 + qJD(2) * t104) * t119, -t107 * t132 + (-qJ(2) * t125 + qJD(2) * t106) * t119, t19 * qJD(3) - t90 * t62 + t86 * t68 - t81 * t70 + t80 * t78, -t20 * qJD(3) - t90 * t63 - t81 * t71 + t86 * t72 + t80 * t77, -t19 * t72 - t20 * t68 + t36 * t63 + t37 * t62 - t111, t13 * t36 + t14 * t37 + t17 * t19 + t18 * t20 + t80 * t90 + t81 * t86, t1 * t31 + t112 * t6, -t1 * t30 - t112 * t7 - t31 * t2 - t6 * t26, t143, -t141, 0, t41 * t26 + t56 * t2 + t32 * t30 + t38 * t7 + (-t103 * t12 + t105 * t11 + (-t103 * t21 - t105 * t22) * qJD(5)) * t96, t41 * t112 + t56 * t1 + t32 * t31 + t38 * t6 - (t103 * t11 + t105 * t12 + (-t103 * t22 + t105 * t21) * qJD(5)) * t96; 0, 0, 0, 0, -t109, -t131, 0, 0, 0, 0, 0, t128 * t104, t128 * t106, -qJD(1) * t68 - t71 * qJD(3), -qJD(1) * t72 + t70 * qJD(3), t78 * t62 + t77 * t63 + t70 * t68 + t71 * t72, -t81 * qJD(1) + t111, 0, 0, 0, 0, 0, -qJD(1) * t26 + t143, -qJD(1) * t112 - t141; 0, 0, 0, 0, 0, 0, t104 * t130, t134 * t109, 0, 0, 0, -qJ(2) * t130, t104 * t131, -t23 * qJD(3) - t68 * t118 - t81 * t72 + t13, t24 * qJD(3) - t72 * t118 + t81 * t68 - t14, (t18 + t23) * t72 + (-t17 + t24) * t68 + (t101 * t62 + t102 * t63) * pkin(3), -t17 * t23 - t18 * t24 + (t101 * t14 + t102 * t13 - t81 * t126) * pkin(3), t151, t149, t152, t150, 0, -t44 * t26 - (-t103 * t16 + t105 * t15) * t96 + ((-t103 * t91 - t105 * t139) * t96 + t113) * qJD(5) + t145, -t105 * t5 - t103 * t4 - t44 * t112 + (t103 * t15 + t105 * t16) * t96 + (-(-t103 * t139 + t105 * t91) * t96 - t105 * t9) * qJD(5) + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 + (t72 - t116) * qJD(3), -0.2e1 * t68 * qJD(3), -t68 ^ 2 - t72 ^ 2, t17 * t72 + t18 * t68 + t80, 0, 0, 0, 0, 0, t2 + t138, t1 - t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t149, t152, t150, 0, t146 * t113 + t145, (-t10 * t96 - t4) * t103 + (-t146 * t9 - t5) * t105 + t148;];
tauc_reg = t3;
