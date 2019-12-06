% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:15
% EndTime: 2019-12-05 17:38:18
% DurationCPUTime: 0.57s
% Computational Cost: add. (923->139), mult. (1888->191), div. (0->0), fcn. (1070->4), ass. (0->99)
t49 = qJD(1) * qJ(2) + qJD(3);
t40 = -pkin(6) * qJD(1) + t49;
t59 = sin(qJ(4));
t88 = t59 * qJD(2);
t61 = cos(qJ(4));
t91 = qJD(4) * t61;
t21 = t40 * t91 + (-pkin(7) * t91 + t88) * qJD(1);
t93 = qJD(1) * t61;
t27 = -pkin(7) * t93 + t61 * t40;
t23 = qJD(4) * pkin(4) + t27;
t60 = cos(qJ(5));
t112 = (qJD(5) * t23 + t21) * t60;
t57 = pkin(1) + qJ(3);
t111 = qJD(1) * t57;
t58 = sin(qJ(5));
t32 = t58 * t61 + t60 * t59;
t28 = t32 * qJD(1);
t42 = pkin(4) * t91 + qJD(3);
t35 = t42 * qJD(1);
t54 = t59 ^ 2;
t55 = t61 ^ 2;
t96 = t54 + t55;
t110 = t96 * qJD(2);
t52 = qJD(4) + qJD(5);
t109 = t28 ^ 2;
t81 = t60 * t93;
t94 = qJD(1) * t59;
t82 = t58 * t94;
t30 = t81 - t82;
t108 = t30 ^ 2;
t56 = -pkin(6) + qJ(2);
t107 = pkin(7) - t56;
t98 = t60 * t61;
t71 = t52 * t98;
t90 = qJD(5) * t58;
t92 = qJD(4) * t59;
t17 = -t58 * t92 - t59 * t90 + t71;
t106 = t17 * t52;
t104 = t30 * t28;
t103 = t30 * t52;
t46 = t59 * pkin(4) + t57;
t31 = t46 * qJD(1) - qJD(2);
t102 = t31 * t30;
t63 = qJD(1) ^ 2;
t101 = t54 * t63;
t26 = -pkin(7) * t94 + t59 * t40;
t100 = t58 * t26;
t99 = t60 * t26;
t85 = qJD(1) * qJD(4);
t78 = t59 * t85;
t97 = -qJD(5) * t82 - t58 * t78;
t62 = qJD(4) ^ 2;
t95 = -t62 - t63;
t89 = t31 * qJD(1);
t41 = -qJD(2) + t111;
t87 = qJD(2) - t41;
t86 = qJD(1) * qJD(2);
t84 = t61 * t63 * t59;
t83 = pkin(4) * t93;
t51 = 0.2e1 * t86;
t80 = 0.2e1 * qJD(3) * qJD(1);
t37 = t107 * t61;
t79 = -pkin(4) * t52 - t23;
t77 = t61 * t85;
t47 = t61 * t86;
t20 = t47 + (pkin(7) * qJD(1) - t40) * t92;
t76 = t60 * t20 - t58 * t21;
t75 = t58 * t20 - t26 * t90;
t74 = t41 + t111;
t72 = t59 * t77;
t16 = t52 * t32;
t12 = t16 * qJD(1);
t33 = -t58 * t59 + t98;
t70 = -t33 * t12 - t30 * t16;
t13 = qJD(1) * t71 + t97;
t69 = t32 * t13 + t17 * t28;
t9 = t58 * t23 + t99;
t36 = t107 * t59;
t19 = -t60 * t36 - t58 * t37;
t18 = t58 * t36 - t60 * t37;
t68 = qJD(2) + t74;
t67 = -t56 * t62 + t80;
t66 = t31 * t28 - t75;
t65 = -t52 * t81 - t97;
t1 = t75 + t112;
t2 = -qJD(5) * t9 + t76;
t8 = t60 * t23 - t100;
t64 = t1 * t32 - t8 * t16 + t9 * t17 + t2 * t33;
t50 = t55 * t63;
t25 = -qJD(4) * t37 + t88;
t24 = t61 * qJD(2) + t107 * t92;
t14 = t16 * t52;
t11 = t60 * t27 - t100;
t10 = -t58 * t27 - t99;
t7 = t108 - t109;
t6 = t65 + t103;
t4 = -t19 * qJD(5) + t60 * t24 - t58 * t25;
t3 = t18 * qJD(5) + t58 * t24 + t60 * t25;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, qJ(2) * t51, 0, 0, 0, 0, 0, 0, 0, t51, t80, t49 * qJD(2) + t41 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t57) * qJD(1), -0.2e1 * t72, 0.2e1 * (t54 - t55) * t85, -t62 * t59, 0.2e1 * t72, -t62 * t61, 0, t67 * t59 + t68 * t91, t67 * t61 - t68 * t92, -t96 * t51, t74 * qJD(3) + (qJD(1) * t56 + t40) * t110, t70, t12 * t32 - t33 * t13 + t16 * t28 - t30 * t17, -t14, t69, -t106, 0, t46 * t13 + t31 * t17 + t42 * t28 + t35 * t32 + t4 * t52, -t46 * t12 - t31 * t16 - t3 * t52 + t42 * t30 + t35 * t33, t18 * t12 - t19 * t13 - t3 * t28 - t4 * t30 - t64, t1 * t19 + t2 * t18 + t9 * t3 + t31 * t42 + t35 * t46 + t8 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t63 * qJ(2), 0, 0, 0, 0, 0, 0, 0, -t63, 0, (-qJD(3) - t49) * qJD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t77, 0.2e1 * t78, t50 + t101, (-t96 * t40 - qJD(3)) * qJD(1), 0, 0, 0, 0, 0, 0, t65 - t103, t28 * t52 + t12, t108 + t109, -t9 * t28 - t8 * t30 - t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t87 * qJD(1), 0, 0, 0, 0, 0, 0, t95 * t59, t95 * t61, 0, (-t41 + t110) * qJD(1), 0, 0, 0, 0, 0, 0, -qJD(1) * t28 - t14, -qJD(1) * t30 - t106, -t69 - t70, t64 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t50 - t101, 0, -t84, 0, 0, -t41 * t93 + t47, -t87 * t94, 0, 0, t104, t7, 0, -t104, t6, 0, -t28 * t83 - t10 * t52 - t102 + (t79 * t58 - t99) * qJD(5) + t76, -t30 * t83 + t11 * t52 + (t79 * qJD(5) - t21) * t60 + t66, (t10 + t9) * t30 + (t11 - t8) * t28 + (t12 * t60 - t13 * t58 + (-t28 * t60 + t30 * t58) * qJD(5)) * pkin(4), -t8 * t10 - t9 * t11 + (-t61 * t89 + t1 * t58 + t2 * t60 + (-t58 * t8 + t60 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t7, 0, -t104, t6, 0, t9 * t52 - t102 + t2, t8 * t52 - t112 + t66, 0, 0;];
tauc_reg = t5;
