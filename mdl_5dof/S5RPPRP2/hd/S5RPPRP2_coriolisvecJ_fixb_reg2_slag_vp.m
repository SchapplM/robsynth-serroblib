% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:24
% EndTime: 2019-12-31 17:49:26
% DurationCPUTime: 0.51s
% Computational Cost: add. (989->129), mult. (2437->158), div. (0->0), fcn. (1710->6), ass. (0->84)
t103 = cos(qJ(4));
t63 = sin(pkin(8));
t65 = cos(pkin(8));
t67 = sin(qJ(4));
t50 = t103 * t63 + t67 * t65;
t109 = t50 * qJD(1);
t87 = t103 * t65;
t82 = qJD(1) * t87;
t54 = qJD(4) * t82;
t91 = qJD(1) * t63;
t88 = t67 * t91;
t31 = qJD(4) * t88 - t54;
t45 = t50 * qJD(4);
t98 = t67 * t63;
t76 = t87 - t98;
t112 = t109 * t45 + t31 * t76;
t32 = qJD(1) * t45;
t40 = -t82 + t88;
t84 = qJD(4) * t103;
t44 = qJD(4) * t98 - t65 * t84;
t96 = -t50 * t32 + t44 * t40;
t114 = t112 - t96;
t113 = t112 + t96;
t108 = t109 ^ 2;
t37 = t40 ^ 2;
t111 = -t37 - t108;
t110 = -t37 + t108;
t57 = sin(pkin(7)) * pkin(1) + qJ(3);
t52 = t57 * qJD(1);
t60 = t65 * qJD(2);
t93 = pkin(6) * qJD(1);
t27 = t60 + (-t52 - t93) * t63;
t30 = t63 * qJD(2) + t65 * t52;
t28 = t65 * t93 + t30;
t10 = t103 * t28 + t67 * t27;
t73 = t50 * qJD(3);
t71 = qJD(1) * t73;
t3 = t10 * qJD(4) + t71;
t105 = pkin(6) + t57;
t46 = t105 * t63;
t47 = t105 * t65;
t77 = -t103 * t46 - t67 * t47;
t107 = t3 * t77;
t106 = t3 * t76;
t7 = qJD(4) * qJ(5) + t10;
t104 = -t10 + t7;
t51 = -cos(pkin(7)) * pkin(1) - pkin(2) - t65 * pkin(3);
t36 = qJD(1) * t51 + qJD(3);
t14 = t40 * pkin(4) - qJ(5) * t109 + t36;
t102 = t14 * t109;
t101 = t109 * t40;
t99 = t67 * t28;
t95 = qJD(3) * t82 + t27 * t84;
t94 = t63 ^ 2 + t65 ^ 2;
t9 = t103 * t27 - t99;
t92 = qJD(5) - t9;
t35 = qJD(4) * t45;
t12 = qJD(3) * t76 + qJD(4) * t77;
t90 = t12 * qJD(4);
t18 = t103 * t47 - t67 * t46;
t13 = qJD(4) * t18 + t73;
t89 = t13 * qJD(4);
t86 = t9 + t99;
t85 = qJD(3) * t91;
t83 = qJD(1) * t94;
t81 = t32 * pkin(4) + t31 * qJ(5);
t80 = (-t63 * t52 + t60) * t63 - t30 * t65;
t78 = -t32 * t76 + t40 * t45;
t75 = t67 * t85 - t95;
t72 = t109 * t13 - t12 * t40 - t18 * t32 + t3 * t50 + t31 * t77;
t70 = 0.2e1 * t109 * qJD(4);
t69 = qJD(3) * t109;
t34 = t44 * qJD(4);
t21 = pkin(4) * t109 + t40 * qJ(5);
t20 = t54 + (t40 - t88) * qJD(4);
t19 = -t54 + (t40 + t88) * qJD(4);
t16 = -pkin(4) * t76 - t50 * qJ(5) + t51;
t15 = t45 * pkin(4) + t44 * qJ(5) - t50 * qJD(5);
t8 = -qJD(5) * t109 + t81;
t6 = -qJD(4) * pkin(4) + t92;
t5 = -t109 * t44 - t31 * t50;
t2 = (-qJD(4) * t28 - t85) * t67 + t95;
t1 = (qJD(5) - t99) * qJD(4) - t75;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t83, (t57 * t83 - t80) * qJD(3), t5, -t114, -t34, t78, -t35, 0, t51 * t32 + t36 * t45 - t89, -t51 * t31 - t36 * t44 - t90, -t10 * t45 + t2 * t76 + t9 * t44 + t72, t10 * t12 - t9 * t13 + t2 * t18 - t107, t5, -t34, t114, 0, t35, t78, t14 * t45 + t15 * t40 + t16 * t32 - t76 * t8 - t89, t1 * t76 - t6 * t44 - t7 * t45 + t72, -t109 * t15 + t14 * t44 + t16 * t31 - t8 * t50 + t90, t1 * t18 + t7 * t12 + t6 * t13 + t14 * t15 + t8 * t16 - t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t34, t113, -t10 * t44 + t2 * t50 - t9 * t45 - t106, 0, 0, 0, 0, 0, 0, -t35, t113, -t34, t1 * t50 - t7 * t44 + t6 * t45 - t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94 * qJD(1) ^ 2, t80 * qJD(1), 0, 0, 0, 0, 0, 0, t70, -t19, t111, t10 * t40 + t109 * t9, 0, 0, 0, 0, 0, 0, t70, t111, t19, t7 * t40 + (-qJD(5) - t6) * t109 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t110, t20, -t101, 0, 0, -t109 * t36 - t69, qJD(4) * t86 + t36 * t40 + t75, 0, 0, t101, t20, -t110, 0, 0, -t101, -t21 * t40 - t102 - t69, pkin(4) * t31 - t32 * qJ(5) + t104 * t109 + (t6 - t92) * t40, -t14 * t40 + t21 * t109 + (0.2e1 * qJD(5) - t86) * qJD(4) - t75, -t3 * pkin(4) + t1 * qJ(5) - t6 * t10 - t14 * t21 + t7 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t20, -qJD(4) ^ 2 - t108, -qJD(4) * t104 + t102 + t71;];
tauc_reg = t4;
