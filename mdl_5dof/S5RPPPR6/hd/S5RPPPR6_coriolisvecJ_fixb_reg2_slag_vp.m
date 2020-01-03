% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:51
% EndTime: 2019-12-31 17:47:54
% DurationCPUTime: 0.84s
% Computational Cost: add. (1154->156), mult. (2866->269), div. (0->0), fcn. (1937->6), ass. (0->102)
t124 = pkin(3) + qJ(2);
t64 = sin(pkin(7));
t105 = qJD(1) * t64;
t67 = sin(qJ(5));
t66 = cos(pkin(7));
t68 = cos(qJ(5));
t109 = t66 * t68;
t63 = sin(pkin(8));
t93 = t63 * t109;
t31 = qJD(1) * t93 - t67 * t105;
t110 = t66 * t67;
t112 = t64 * t68;
t34 = t63 * t110 + t112;
t28 = t34 * qJD(5);
t23 = qJD(1) * t28;
t27 = t34 * qJD(1);
t104 = qJD(1) * t66;
t65 = cos(pkin(8));
t44 = t65 * t104 + qJD(5);
t123 = -t27 * t44 + t23;
t24 = t31 * qJD(5);
t122 = -t31 * t44 + t24;
t98 = t64 * qJD(3);
t43 = -t66 * qJD(4) - t98;
t41 = t43 * qJD(1);
t96 = qJD(1) * qJD(2);
t86 = t64 * t96;
t17 = t65 * t41 + t63 * t86;
t84 = -t64 * qJ(3) - pkin(1);
t37 = (-pkin(2) - qJ(4)) * t66 + t84;
t20 = t37 * qJD(1) + qJD(2);
t97 = qJ(2) * qJD(1);
t49 = t64 * t97 + qJD(3);
t38 = pkin(3) * t105 + t49;
t12 = t65 * t20 + t63 * t38;
t10 = pkin(6) * t105 + t12;
t39 = pkin(3) * t104 + t66 * t97 + qJD(4);
t75 = (pkin(4) * t65 + pkin(6) * t63) * t66;
t15 = qJD(1) * t75 + t39;
t80 = t67 * t10 - t68 * t15;
t85 = t66 * t96;
t1 = -t80 * qJD(5) + t68 * t17 + t67 * t85;
t16 = t63 * t41 - t65 * t86;
t121 = t16 * t63;
t120 = t27 * t31;
t117 = t44 * t65;
t60 = t64 ^ 2;
t69 = qJD(1) ^ 2;
t116 = t60 * t69;
t115 = t63 * t69;
t114 = t64 * t66;
t113 = t64 * t67;
t111 = t65 * t66;
t45 = t124 * t64;
t108 = t65 * t37 + t63 * t45;
t107 = t124 * t66;
t62 = t66 ^ 2;
t106 = t60 + t62;
t103 = qJD(2) * t62;
t102 = qJD(2) * t64;
t101 = qJD(2) * t66;
t100 = qJD(5) * t67;
t99 = qJD(5) * t68;
t95 = qJD(1) * qJD(3);
t94 = t69 * t114;
t92 = 0.2e1 * t103;
t91 = t65 * t105;
t90 = t44 * t100;
t89 = t44 * t99;
t87 = (-t63 ^ 2 - t65 ^ 2) * t69;
t42 = t106 * t69;
t83 = qJ(2) * t96;
t6 = t68 * t10 + t67 * t15;
t81 = t6 * t68 + t67 * t80;
t11 = -t63 * t20 + t65 * t38;
t79 = -t11 * t63 + t12 * t65;
t14 = t64 * pkin(6) + t108;
t19 = t75 + t107;
t8 = t68 * t14 + t67 * t19;
t7 = -t67 * t14 + t68 * t19;
t78 = -t17 * t65 - t121;
t77 = -t63 * t37 + t65 * t45;
t76 = -t66 * pkin(2) + t84;
t35 = t93 - t113;
t73 = qJD(1) * (-t67 * t117 - t27 * t63);
t72 = qJD(1) * (-t68 * t117 - t31 * t63);
t2 = -t6 * qJD(5) - t67 * t17 + t68 * t85;
t71 = t23 * t67 + t24 * t68 + (-t27 * t67 - t31 * t68) * qJD(5);
t9 = -pkin(4) * t105 - t11;
t70 = t9 * t105 + t1 * t68 - t2 * t67 + (-t6 * t67 + t68 * t80) * qJD(5);
t51 = t60 * t83;
t36 = 0.2e1 * t106 * t96;
t33 = t76 * qJD(1) + qJD(2);
t30 = (t63 * t112 + t110) * qJD(1);
t29 = t35 * qJD(5);
t26 = (-t63 * t113 + t109) * qJD(1);
t22 = t63 * t102 + t65 * t43;
t21 = -t65 * t102 + t63 * t43;
t13 = -t64 * pkin(4) - t77;
t4 = -t8 * qJD(5) + t68 * t101 - t67 * t22;
t3 = t7 * qJD(5) + t67 * t101 + t68 * t22;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0.2e1 * t62 * t83 + 0.2e1 * t51, 0, 0, 0, 0, 0, 0, t36, -0.2e1 * t95 * t114, 0.2e1 * t60 * t95, t51 + (t49 * qJD(2) - t33 * qJD(3)) * t64 + (qJ(2) * t92 - t76 * t98) * qJD(1), 0, 0, 0, 0, 0, 0, -t16 * t64 + (-t21 * t64 + t65 * t92) * qJD(1), -t17 * t64 + (-0.2e1 * t63 * t103 - t22 * t64) * qJD(1), ((-t21 * t63 - t22 * t65) * qJD(1) + t78) * t66, t17 * t108 + t12 * t22 - t16 * t77 - t11 * t21 + (qJD(1) * t107 + t39) * t101, -t23 * t35 - t31 * t28, t23 * t34 - t35 * t24 + t28 * t27 - t31 * t29, t23 * t111 + t28 * t44, t24 * t34 + t27 * t29, t24 * t111 + t29 * t44, 0, t2 * t111 - t13 * t24 - t16 * t34 - t21 * t27 - t9 * t29 + t4 * t44, -t1 * t111 + t13 * t23 - t16 * t35 - t21 * t31 + t9 * t28 - t3 * t44, t1 * t34 + t2 * t35 - t7 * t23 + t8 * t24 + t3 * t27 + t28 * t80 + t6 * t29 + t4 * t31, t1 * t8 + t16 * t13 + t2 * t7 + t9 * t21 + t6 * t3 - t4 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -qJ(2) * t42, 0, 0, 0, 0, 0, 0, -t42, 0, 0, -t62 * t69 * qJ(2) + (-qJD(3) - t49) * t105, 0, 0, 0, 0, 0, 0, -t65 * t42, t106 * t115, 0, (-t39 * t66 + (-t11 * t65 - t12 * t63) * t64) * qJD(1) - t78, 0, 0, 0, 0, 0, 0, -t27 * t91 - t63 * t24 + (-t65 * t99 - t26) * t44, -t31 * t91 + t63 * t23 + (t65 * t100 + t30) * t44, -t26 * t31 - t30 * t27 + t71 * t65, t26 * t80 - t6 * t30 + t70 * t65 + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t116, (qJD(2) + t33) * t105, 0, 0, 0, 0, 0, 0, -t60 * t115, -t65 * t116, t87 * t114, t79 * t105 - t16 * t65 + t17 * t63, 0, 0, 0, 0, 0, 0, t65 * t24 - t63 * t89 + t64 * t73, -t65 * t23 + t63 * t90 + t64 * t72, (t27 * t68 - t31 * t67) * t91 + t71 * t63, (t81 * t105 - t16) * t65 + t70 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63 * t94, -t65 * t94, t62 * t87, (qJD(2) + t79) * t104, 0, 0, 0, 0, 0, 0, t66 * t73 - t90, t66 * t72 - t89, t122 * t67 - t123 * t68, t1 * t67 + t2 * t68 + t81 * qJD(5) + (t63 * t9 + t81 * t65) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t27 ^ 2 + t31 ^ 2, t123, -t120, t122, 0, t9 * t31 + t6 * t44 + t2, -t9 * t27 - t44 * t80 - t1, 0, 0;];
tauc_reg = t5;
