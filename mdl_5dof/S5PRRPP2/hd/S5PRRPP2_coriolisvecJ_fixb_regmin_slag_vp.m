% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:32
% EndTime: 2021-01-15 15:32:37
% DurationCPUTime: 0.67s
% Computational Cost: add. (1026->164), mult. (2653->232), div. (0->0), fcn. (1792->6), ass. (0->100)
t77 = sin(qJ(3));
t79 = cos(qJ(3));
t124 = qJ(4) + pkin(6);
t99 = qJD(3) * t124;
t44 = t79 * qJD(4) - t77 * t99;
t75 = sin(pkin(8));
t76 = cos(pkin(8));
t88 = -t77 * qJD(4) - t79 * t99;
t127 = t76 * t79;
t55 = t75 * t77 - t127;
t80 = cos(qJ(2));
t91 = t55 * t80;
t122 = qJD(1) * t91 + t76 * t44 + t75 * t88;
t109 = (qJD(2) * qJD(3));
t138 = -2 * t109;
t56 = t75 * t79 + t76 * t77;
t50 = t56 * qJD(2);
t46 = t50 ^ 2;
t104 = qJD(2) * t127;
t118 = qJD(2) * t77;
t47 = t75 * t118 - t104;
t137 = -t47 ^ 2 - t46;
t78 = sin(qJ(2));
t81 = qJD(3) ^ 2;
t82 = qJD(2) ^ 2;
t136 = (t81 + t82) * t78;
t49 = t56 * qJD(3);
t114 = qJD(3) * t79;
t115 = qJD(3) * t77;
t52 = t76 * t114 - t75 * t115;
t113 = t78 * qJD(1);
t135 = pkin(3) * t115 - t113;
t134 = pkin(3) * t77;
t110 = qJ(4) * qJD(3);
t63 = qJD(2) * pkin(6) + t113;
t112 = t80 * qJD(1);
t95 = qJD(4) + t112;
t26 = -t63 * t115 + (-t77 * t110 + t95 * t79) * qJD(2);
t83 = -t63 * t114 + (-t79 * t110 - t95 * t77) * qJD(2);
t2 = t75 * t26 - t76 * t83;
t103 = t124 * t77;
t61 = t124 * t79;
t28 = t76 * t103 + t75 * t61;
t133 = t2 * t28;
t40 = t56 * t78;
t132 = t2 * t40;
t72 = -t79 * pkin(3) - pkin(2);
t53 = t72 * qJD(2) + qJD(4) - t112;
t9 = t47 * pkin(4) - t50 * qJ(5) + t53;
t131 = t9 * t50;
t97 = qJ(4) * qJD(2) + t63;
t43 = t97 * t79;
t130 = t75 * t43;
t31 = t76 * t43;
t126 = t81 * t77;
t125 = t81 * t79;
t123 = -t56 * t112 + t75 * t44 - t76 * t88;
t3 = t76 * t26 + t75 * t83;
t42 = t97 * t77;
t33 = qJD(3) * pkin(3) - t42;
t12 = t75 * t33 + t31;
t102 = t77 * t109;
t54 = pkin(3) * t102 + qJD(2) * t113;
t121 = t77 ^ 2 - t79 ^ 2;
t119 = qJD(2) * pkin(2);
t117 = qJD(2) * t78;
t15 = -t76 * t42 - t130;
t111 = qJD(5) - t15;
t107 = pkin(3) * t118;
t101 = t79 * t109;
t100 = -t49 * pkin(4) + t52 * qJ(5) + t56 * qJD(5) - t135;
t96 = t80 * t138;
t11 = t76 * t33 - t130;
t94 = qJD(2) * t119;
t14 = -t75 * t42 + t31;
t93 = t14 * qJD(3) - t2;
t16 = -t80 * t50 - t52 * t78;
t38 = qJD(2) * t49;
t92 = t16 * qJD(3) + t47 * t117 - t80 * t38;
t17 = -qJD(2) * t91 - t78 * t49;
t62 = t75 * t102;
t39 = t76 * t101 - t62;
t41 = t55 * t78;
t90 = -t16 * t50 - t17 * t47 + t41 * t38 + t40 * t39;
t89 = t38 * pkin(4) - t39 * qJ(5) + t54;
t87 = -0.2e1 * qJD(3) * t119;
t86 = t17 * qJD(3) - t50 * t117 + t80 * t39;
t85 = 0.2e1 * t50 * qJD(3);
t29 = -t75 * t103 + t76 * t61;
t84 = -t122 * t47 + t123 * t50 + t2 * t56 + t28 * t39 - t29 * t38;
t70 = -t76 * pkin(3) - pkin(4);
t68 = t75 * pkin(3) + qJ(5);
t27 = t55 * pkin(4) - t56 * qJ(5) + t72;
t25 = -t62 + (-t47 + t104) * qJD(3);
t18 = t50 * pkin(4) + t47 * qJ(5) + t107;
t10 = qJD(3) * qJ(5) + t12;
t8 = -qJD(3) * pkin(4) + qJD(5) - t11;
t4 = -t50 * qJD(5) + t89;
t1 = qJD(3) * qJD(5) + t3;
t5 = [0, 0, -t82 * t78, -t82 * t80, 0, 0, 0, 0, 0, -t79 * t136 + t77 * t96, t77 * t136 + t79 * t96, t92, -t86, t90, t11 * t16 + t53 * t117 + t12 * t17 - t3 * t41 - t54 * t80 + t132, t92, t90, t86, -t1 * t41 + t10 * t17 + t9 * t117 - t8 * t16 - t4 * t80 + t132; 0, 0, 0, 0, 0.2e1 * t77 * t101, t121 * t138, t125, -t126, 0, -pkin(6) * t125 + t77 * t87, pkin(6) * t126 + t79 * t87, -t47 * t113 + t72 * t38 + t53 * t49 + t54 * t55 + (t47 * t134 - t123) * qJD(3), -t50 * t113 + t72 * t39 + t53 * t52 + t54 * t56 + (t50 * t134 - t122) * qJD(3), -t11 * t52 - t12 * t49 - t3 * t55 + t84, -t123 * t11 + t122 * t12 + t135 * t53 + t3 * t29 + t54 * t72 + t133, -t123 * qJD(3) - t100 * t47 + t27 * t38 + t4 * t55 + t9 * t49, -t1 * t55 - t10 * t49 + t8 * t52 + t84, t122 * qJD(3) + t100 * t50 - t27 * t39 - t4 * t56 - t9 * t52, t1 * t29 + t122 * t10 - t100 * t9 + t123 * t8 + t4 * t27 + t133; 0, 0, 0, 0, -t77 * t82 * t79, t121 * t82, 0, 0, 0, t77 * t94, t79 * t94, -t47 * t107 - t53 * t50 + t93, t15 * qJD(3) - t50 * t107 + t53 * t47 - t3, (t12 - t14) * t50 + (-t11 + t15) * t47 + (-t38 * t75 - t39 * t76) * pkin(3), t11 * t14 - t12 * t15 + (-t53 * t118 - t2 * t76 + t3 * t75) * pkin(3), -t18 * t47 - t131 + t93, -t68 * t38 + t70 * t39 + (t10 - t14) * t50 + (t8 - t111) * t47, t18 * t50 - t9 * t47 + (0.2e1 * qJD(5) - t15) * qJD(3) + t3, t1 * t68 + t111 * t10 - t8 * t14 - t9 * t18 + t2 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t25, t137, t11 * t50 + t12 * t47 + t54, t85, t137, -t25, t10 * t47 + (-qJD(5) - t8) * t50 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t47, -t62 + (t47 + t104) * qJD(3), -t46 - t81, -t10 * qJD(3) + t131 + t2;];
tauc_reg = t5;
