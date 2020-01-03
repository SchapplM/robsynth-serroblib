% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:17
% EndTime: 2019-12-31 17:55:19
% DurationCPUTime: 0.65s
% Computational Cost: add. (1021->137), mult. (2231->161), div. (0->0), fcn. (1462->4), ass. (0->90)
t61 = sin(pkin(7));
t64 = sin(qJ(4));
t62 = cos(pkin(7));
t65 = cos(qJ(4));
t98 = t65 * t62;
t42 = -t64 * t61 + t98;
t119 = t42 * qJD(3);
t86 = qJD(1) * t98;
t94 = qJD(1) * t61;
t88 = t64 * t94;
t36 = t86 - t88;
t60 = qJD(1) * qJ(2);
t55 = qJD(3) + t60;
t45 = pkin(3) * t94 + t55;
t41 = t61 * t65 + t62 * t64;
t74 = qJD(1) * t41;
t13 = pkin(4) * t74 - qJ(5) * t36 + t45;
t77 = t119 * qJD(1);
t118 = t13 * t36 + t77;
t38 = t41 * qJD(4);
t29 = t38 * qJD(4);
t117 = -qJD(1) * t74 - t29;
t116 = 0.2e1 * t74;
t47 = qJD(4) * t88;
t26 = qJD(4) * t86 - t47;
t92 = qJD(4) * t65;
t93 = qJD(4) * t64;
t39 = -t61 * t93 + t62 * t92;
t78 = t26 * t41 + t39 * t74;
t109 = t36 ^ 2;
t31 = t74 ^ 2;
t115 = -t31 - t109;
t114 = -t31 + t109;
t63 = -pkin(1) - qJ(3);
t113 = t63 * qJD(1);
t96 = t61 ^ 2 + t62 ^ 2;
t112 = t96 * qJD(3);
t106 = -pkin(6) + t63;
t43 = t106 * t61;
t44 = t106 * t62;
t22 = t43 * t64 - t44 * t65;
t75 = t41 * qJD(3);
t10 = -qJD(4) * t22 - t75;
t23 = t43 * t65 + t44 * t64;
t11 = qJD(4) * t23 + t119;
t25 = qJD(1) * t38;
t111 = -t10 * t74 + t11 * t36 - t22 * t25 - t23 * t26;
t110 = qJD(4) * (-t36 + t86) - t47;
t59 = qJD(1) * qJD(2);
t84 = 0.2e1 * t59;
t48 = qJD(2) + t113;
t82 = -pkin(6) * qJD(1) + t48;
t27 = t82 * t61;
t28 = t82 * t62;
t15 = t27 * t65 + t28 * t64;
t4 = t15 * qJD(4) + t77;
t108 = t4 * t22;
t107 = t4 * t42;
t101 = t36 * t74;
t100 = t64 * t27;
t12 = qJD(4) * qJ(5) + t15;
t97 = t12 - t15;
t52 = pkin(3) * t61 + qJ(2);
t91 = t10 * qJD(4);
t90 = t11 * qJD(4);
t30 = t39 * qJD(4);
t14 = t28 * t65 - t100;
t89 = qJD(5) - t14;
t83 = t14 + t100;
t81 = qJD(1) * t96;
t80 = qJD(1) * t36 + t30;
t2 = -t25 * t42 - t36 * t38;
t76 = t26 * pkin(4) + t25 * qJ(5) + t59;
t73 = -t2 - t78;
t72 = -t47 + (t36 + t86) * qJD(4);
t24 = t28 * t92;
t70 = -qJD(1) * t75 + t24;
t1 = (qJD(5) - t100) * qJD(4) + t70;
t9 = -qJD(4) * pkin(4) + t89;
t71 = t1 * t41 + t12 * t39 + t38 * t9 - t107;
t3 = -t27 * t93 + t70;
t69 = -t14 * t38 + t15 * t39 + t3 * t41 - t107;
t68 = t25 * t41 - t26 * t42 - t36 * t39 + t38 * t74;
t66 = qJD(1) ^ 2;
t19 = pkin(4) * t36 + qJ(5) * t74;
t18 = pkin(4) * t41 - qJ(5) * t42 + t52;
t16 = t116 * qJD(4);
t7 = pkin(4) * t39 + qJ(5) * t38 - qJD(5) * t42 + qJD(2);
t6 = -qJD(5) * t36 + t76;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, qJ(2) * t84, 0, 0, 0, 0, 0, 0, t61 * t84, t62 * t84, 0.2e1 * qJD(3) * t81, (t55 + t60) * qJD(2) + (-t48 - t113) * t112, t2, t68, -t29, t78, -t30, 0, qJD(2) * t116 + t52 * t26 + t45 * t39 - t90, -t91 - t52 * t25 - t45 * t38 + (qJD(1) * t42 + t36) * qJD(2), -t69 + t111, t15 * t10 - t14 * t11 + t108 + t3 * t23 + (qJD(1) * t52 + t45) * qJD(2), t2, -t29, -t68, 0, t30, t78, t13 * t39 + t18 * t26 + t41 * t6 + t7 * t74 - t90, -t71 + t111, t13 * t38 + t18 * t25 - t36 * t7 - t42 * t6 + t91, t1 * t23 + t10 * t12 + t11 * t9 + t13 * t7 + t18 * t6 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t66 * qJ(2), 0, 0, 0, 0, 0, 0, -t66 * t61, -t66 * t62, 0, (-t55 - t112) * qJD(1), 0, 0, 0, 0, 0, 0, t117, -t80, t73, -qJD(1) * t45 + t69, 0, 0, 0, 0, 0, 0, t117, t73, t80, -qJD(1) * t13 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96 * t66, t48 * t81 + t59, 0, 0, 0, 0, 0, 0, t72, -t16, t115, t14 * t36 + t15 * t74 + t59, 0, 0, 0, 0, 0, 0, t72, t115, t16, t12 * t74 + (-qJD(5) - t9) * t36 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t114, 0, -t101, -t110, 0, -t45 * t36 - t77, qJD(4) * t83 - t24 + (qJD(3) + t45) * t74, 0, 0, t101, 0, -t114, 0, t110, -t101, -t19 * t74 - t118, pkin(4) * t25 - t26 * qJ(5) + t97 * t36 + (t9 - t89) * t74, -t13 * t74 + t19 * t36 + (0.2e1 * qJD(5) - t83) * qJD(4) + t70, -t4 * pkin(4) + t1 * qJ(5) + t12 * t89 - t13 * t19 - t9 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, -qJD(4) ^ 2 - t109, -qJD(4) * t97 + t118;];
tauc_reg = t5;
