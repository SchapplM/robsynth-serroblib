% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:42
% EndTime: 2019-12-31 17:42:44
% DurationCPUTime: 0.58s
% Computational Cost: add. (1388->119), mult. (3060->178), div. (0->0), fcn. (2247->8), ass. (0->96)
t61 = sin(pkin(9));
t64 = sin(qJ(3));
t105 = t61 * t64;
t65 = sin(qJ(2));
t67 = cos(qJ(3));
t103 = t67 * t65;
t68 = cos(qJ(2));
t44 = t64 * t68 + t103;
t41 = t44 * qJD(1);
t43 = -t64 * t65 + t67 * t68;
t42 = t43 * qJD(1);
t62 = cos(pkin(9));
t95 = pkin(2) * qJD(3);
t100 = t61 * t41 - t62 * t42 + (t62 * t67 - t105) * t95;
t58 = qJD(2) + qJD(3);
t118 = t100 * t58;
t91 = t68 * qJD(1);
t51 = qJD(2) * pkin(2) + t91;
t94 = qJD(1) * t65;
t89 = t64 * t94;
t33 = t67 * t51 - t89;
t30 = t58 * pkin(3) + t33;
t34 = t64 * t51 + t67 * t94;
t31 = t62 * t34;
t16 = t61 * t30 + t31;
t14 = t58 * pkin(7) + t16;
t63 = sin(qJ(5));
t66 = cos(qJ(5));
t10 = t63 * qJD(4) + t66 * t14;
t86 = qJD(2) * t91;
t93 = qJD(3) * t51;
t99 = t58 * t89;
t20 = t67 * (t86 + t93) - t99;
t76 = t44 * qJD(2);
t71 = (-qJD(3) * t103 - t76) * qJD(1);
t21 = -t64 * t93 + t71;
t6 = t62 * t20 + t61 * t21;
t9 = t66 * qJD(4) - t63 * t14;
t2 = qJD(5) * t9 + t66 * t6;
t3 = -qJD(5) * t10 - t63 * t6;
t117 = t2 * t66 - t3 * t63 + (-t10 * t63 - t66 * t9) * qJD(5);
t5 = t61 * t20 - t62 * t21;
t104 = t62 * t64;
t101 = t62 * t41 + t61 * t42 - (t61 * t67 + t104) * t95;
t87 = t101 * t58;
t116 = t87 - t5;
t84 = t10 * t66 - t63 * t9;
t25 = -t62 * t43 + t61 * t44;
t115 = t5 * t25;
t27 = t58 * t43;
t28 = -t44 * qJD(3) - t76;
t7 = t61 * t27 - t62 * t28;
t113 = t7 * t58;
t8 = t62 * t27 + t61 * t28;
t112 = t8 * t58;
t106 = t61 * t34;
t15 = t62 * t30 - t106;
t13 = -t58 * pkin(4) - t15;
t92 = qJD(5) * t66;
t111 = t13 * t92 + t5 * t63;
t17 = t61 * t33 + t31;
t109 = t17 * t58;
t18 = t62 * t33 - t106;
t108 = t18 * t58;
t55 = t67 * pkin(2) + pkin(3);
t98 = pkin(2) * t104 + t61 * t55;
t38 = pkin(7) + t98;
t69 = qJD(5) ^ 2;
t107 = t38 * t69;
t102 = t69 * t63;
t59 = t63 ^ 2;
t60 = t66 ^ 2;
t97 = t59 - t60;
t96 = t59 + t60;
t57 = t58 ^ 2;
t90 = t63 * t57 * t66;
t88 = -t13 * t58 - t6;
t85 = t63 * t58 * t92;
t26 = t61 * t43 + t62 * t44;
t83 = t26 * t69 + t113;
t53 = t61 * pkin(3) + pkin(7);
t82 = t53 * t69 - t109;
t81 = qJD(5) * (t25 * t58 - t8);
t80 = (-pkin(2) * t58 - t51) * qJD(3);
t54 = -t62 * pkin(3) - pkin(4);
t79 = qJD(5) * (t54 * t58 + t18);
t78 = -pkin(2) * t105 + t62 * t55;
t37 = -pkin(4) - t78;
t77 = qJD(5) * (t37 * t58 - t100);
t70 = qJD(2) ^ 2;
t56 = t69 * t66;
t47 = -0.2e1 * t85;
t46 = 0.2e1 * t85;
t35 = -0.2e1 * t97 * t58 * qJD(5);
t11 = t13 * qJD(5) * t63;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70 * t65, -t70 * t68, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t58, -t27 * t58, 0, t20 * t44 + t21 * t43 + t34 * t27 + t33 * t28, 0, 0, 0, 0, 0, 0, -t113, -t112, 0, -t15 * t7 + t16 * t8 + t6 * t26 + t115, 0, 0, 0, 0, 0, 0, t63 * t81 - t66 * t83, t63 * t83 + t66 * t81, t96 * t112, t117 * t26 + t13 * t7 + t8 * t84 + t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t58 + t64 * t80 + t71, t42 * t58 + (t80 - t86) * t67 + t99, 0, t33 * t41 - t34 * t42 + (t20 * t64 + t21 * t67 + (-t33 * t64 + t34 * t67) * qJD(3)) * pkin(2), 0, 0, 0, 0, 0, 0, t116, -t6 - t118, 0, t100 * t16 + t101 * t15 - t5 * t78 + t6 * t98, t46, t35, t56, t47, -t102, 0, t11 + t63 * t77 + (-t107 + t116) * t66, (t107 - t87) * t63 + t66 * t77 + t111, t96 * t118 + t117, t84 * t100 - t101 * t13 + t117 * t38 + t5 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t58 + t21, t33 * t58 - t20, 0, 0, 0, 0, 0, 0, 0, 0, -t5 + t109, -t6 + t108, 0, t15 * t17 - t16 * t18 + (-t5 * t62 + t6 * t61) * pkin(3), t46, t35, t56, t47, -t102, 0, t11 + t63 * t79 + (-t5 - t82) * t66, t63 * t82 + t66 * t79 + t111, -t96 * t108 + t117, t117 * t53 - t13 * t17 - t18 * t84 + t5 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t56, 0, qJD(5) * t84 + t2 * t63 + t3 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, t97 * t57, 0, t90, 0, 0, t88 * t63, t88 * t66, 0, 0;];
tauc_reg = t1;
