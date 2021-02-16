% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:23
% EndTime: 2021-01-15 14:30:26
% DurationCPUTime: 0.62s
% Computational Cost: add. (852->139), mult. (2305->194), div. (0->0), fcn. (1496->4), ass. (0->100)
t65 = sin(qJ(3));
t66 = sin(qJ(2));
t67 = cos(qJ(3));
t68 = cos(qJ(2));
t44 = t65 * t68 + t67 * t66;
t102 = qJD(1) * t44;
t103 = t102 * qJ(4);
t122 = pkin(5) + pkin(6);
t51 = t122 * t68;
t47 = qJD(1) * t51;
t33 = t65 * t47;
t107 = qJD(2) * pkin(2);
t50 = t122 * t66;
t45 = qJD(1) * t50;
t39 = -t45 + t107;
t90 = t67 * t39 - t33;
t8 = -t103 + t90;
t62 = qJD(2) + qJD(3);
t20 = t62 * t44;
t18 = t20 * qJD(1);
t100 = qJD(3) * t65;
t92 = qJD(2) * t122;
t83 = qJD(1) * t92;
t40 = t66 * t83;
t41 = t68 * t83;
t77 = -(qJD(3) * t39 - t40) * t67 + t47 * t100 + t65 * t41;
t125 = -t18 * qJ(4) - t77;
t97 = qJD(1) * qJD(2);
t124 = -0.2e1 * t97;
t123 = t102 ^ 2;
t5 = t62 * pkin(3) + t8;
t121 = t5 - t8;
t120 = pkin(2) * t62;
t114 = t67 * t68;
t93 = qJD(1) * t114;
t101 = qJD(1) * t66;
t94 = t65 * t101;
t29 = -t93 + t94;
t119 = t29 * t62;
t118 = t102 * t29;
t61 = -t68 * pkin(2) - pkin(1);
t49 = t61 * qJD(1);
t116 = t49 * t102;
t115 = t65 * t66;
t37 = t67 * t47;
t70 = qJD(1) ^ 2;
t113 = t68 * t70;
t69 = qJD(2) ^ 2;
t112 = t69 * t66;
t111 = t69 * t68;
t110 = -t67 * t45 - t33;
t91 = t68 * t97;
t109 = -qJD(3) * t93 - t67 * t91;
t108 = t66 ^ 2 - t68 ^ 2;
t82 = t62 * t115;
t17 = qJD(1) * t82 + t109;
t106 = t17 * qJ(4);
t104 = t29 * qJ(4);
t99 = qJD(3) * t67;
t86 = t29 * pkin(3) + qJD(4);
t21 = t49 + t86;
t98 = qJD(4) + t21;
t96 = t66 * t107;
t95 = pkin(2) * t101;
t13 = t18 * pkin(3) + qJD(2) * t95;
t89 = t65 * t40 - t67 * t41;
t87 = t65 * t45 - t37;
t84 = pkin(1) * t124;
t81 = -t65 * t39 - t37;
t80 = t65 * t50 - t67 * t51;
t46 = t66 * t92;
t48 = t68 * t92;
t79 = -t51 * t100 - t67 * t46 - t65 * t48 - t50 * t99;
t78 = -t62 * t94 - t109;
t76 = t81 * qJD(3) + t89;
t75 = t80 * qJD(3) + t65 * t46 - t67 * t48;
t74 = t49 * t29 + t77;
t73 = t76 + t106;
t72 = (-t37 + (-t39 - t120) * t65) * qJD(3) + t89;
t71 = t98 * t29 - t125;
t60 = t67 * pkin(2) + pkin(3);
t52 = t99 * t120;
t43 = -t114 + t115;
t28 = t29 ^ 2;
t23 = t43 * pkin(3) + t61;
t22 = pkin(3) * t102 + t95;
t19 = -qJD(2) * t114 - t68 * t99 + t82;
t16 = t20 * pkin(3) + t96;
t15 = -t43 * qJ(4) - t80;
t14 = -t44 * qJ(4) - t67 * t50 - t65 * t51;
t12 = -t28 + t123;
t11 = -t103 + t110;
t10 = t87 + t104;
t9 = -t81 - t104;
t6 = t78 + t119;
t4 = t19 * qJ(4) - t44 * qJD(4) + t75;
t3 = -t20 * qJ(4) - t43 * qJD(4) + t79;
t2 = -qJD(4) * t102 + t73;
t1 = -t29 * qJD(4) + t125;
t7 = [0, 0, 0, 0.2e1 * t66 * t91, t108 * t124, t111, -t112, 0, -pkin(5) * t111 + t66 * t84, pkin(5) * t112 + t68 * t84, -t102 * t19 - t17 * t44, -t102 * t20 + t17 * t43 - t44 * t18 + t19 * t29, -t19 * t62, -t20 * t62, 0, t61 * t18 + t49 * t20 + t75 * t62 + (qJD(1) * t43 + t29) * t96, 0.2e1 * t102 * t96 - t61 * t17 - t49 * t19 - t79 * t62, t13 * t43 + t16 * t29 + t23 * t18 + t21 * t20 + t4 * t62, t102 * t16 + t13 * t44 - t23 * t17 - t21 * t19 - t3 * t62, -t1 * t43 - t102 * t4 + t14 * t17 - t15 * t18 + t5 * t19 - t2 * t44 - t9 * t20 - t3 * t29, t1 * t15 + t13 * t23 + t2 * t14 + t21 * t16 + t9 * t3 + t5 * t4; 0, 0, 0, -t66 * t113, t108 * t70, 0, 0, 0, t70 * pkin(1) * t66, pkin(1) * t113, t118, t12, t6, 0, 0, -t29 * t95 - t87 * t62 - t116 + t72, -t102 * t95 + t110 * t62 - t52 + t74, -t10 * t62 - t102 * t98 - t22 * t29 + t106 + t72, -t102 * t22 + t11 * t62 - t52 + t71, t60 * t17 + (t10 + t9) * t102 + (t11 - t5) * t29 + (-t18 * t65 + (t102 * t65 - t29 * t67) * qJD(3)) * pkin(2), -t5 * t10 - t9 * t11 + t2 * t60 - t21 * t22 + (t1 * t65 + (-t5 * t65 + t67 * t9) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t12, t6, 0, 0, -t81 * t62 - t116 + t76, t90 * t62 + t74, t9 * t62 + (-t21 - t86) * t102 + t73, -t123 * pkin(3) + t8 * t62 + t71, t17 * pkin(3) - t121 * t29, t121 * t9 + (-t102 * t21 + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t102 + t18, t78 - t119, -t28 - t123, t102 * t5 + t9 * t29 + t13;];
tauc_reg = t7;
