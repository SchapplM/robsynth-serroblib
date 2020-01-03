% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:39
% EndTime: 2019-12-31 16:32:41
% DurationCPUTime: 0.43s
% Computational Cost: add. (654->111), mult. (1727->169), div. (0->0), fcn. (1102->4), ass. (0->86)
t55 = sin(qJ(3));
t83 = qJD(3) * pkin(3);
t76 = t55 * t83;
t100 = 0.2e1 * t76;
t79 = (qJD(2) * qJD(3));
t99 = -2 * t79;
t97 = -pkin(6) - pkin(5);
t70 = qJD(3) * t97;
t35 = t55 * t70;
t57 = cos(qJ(3));
t50 = t57 * qJD(1);
t47 = qJD(3) * t50;
t23 = qJD(2) * t35 + t47;
t71 = qJD(2) * t97;
t26 = t55 * t71 + t50;
t25 = t26 + t83;
t56 = cos(qJ(4));
t98 = (qJD(4) * t25 + t23) * t56;
t54 = sin(qJ(4));
t34 = t54 * t57 + t56 * t55;
t30 = qJD(2) * t34;
t51 = qJD(3) + qJD(4);
t18 = t51 * t34;
t15 = t18 * qJD(2);
t91 = t54 * t55;
t63 = t51 * t91;
t89 = t56 * t57;
t17 = -t51 * t89 + t63;
t81 = qJD(2) * t57;
t73 = t56 * t81;
t82 = qJD(2) * t55;
t74 = t54 * t82;
t28 = -t73 + t74;
t96 = -t34 * t15 + t17 * t28;
t95 = t17 * t51;
t94 = t30 * t28;
t49 = -t57 * pkin(3) - pkin(2);
t39 = qJD(2) * t49;
t93 = t39 * t30;
t41 = t97 * t57;
t80 = t55 * qJD(1);
t27 = -qJD(2) * t41 + t80;
t92 = t54 * t27;
t90 = t56 * t27;
t59 = qJD(2) ^ 2;
t88 = t57 * t59;
t58 = qJD(3) ^ 2;
t87 = t58 * t55;
t86 = t58 * t57;
t69 = t57 * t79;
t85 = -qJD(4) * t73 - t56 * t69;
t84 = t55 ^ 2 - t57 ^ 2;
t78 = t55 * t88;
t77 = pkin(3) * t82;
t75 = pkin(5) * t82;
t72 = -pkin(3) * t51 - t25;
t24 = (t57 * t71 - t80) * qJD(3);
t68 = -t54 * t23 + t56 * t24;
t67 = -qJD(4) * t92 + t54 * t24;
t65 = pkin(2) * t99;
t64 = t55 * t69;
t14 = qJD(2) * t63 + t85;
t33 = -t89 + t91;
t62 = -t33 * t14 + t30 * t18;
t11 = t54 * t25 + t90;
t40 = t97 * t55;
t19 = t56 * t40 + t54 * t41;
t20 = t54 * t40 - t56 * t41;
t61 = t39 * t28 - t67;
t38 = pkin(5) * t81 + t80;
t2 = -qJD(4) * t11 + t68;
t31 = -qJD(3) * t75 + t47;
t32 = t38 * qJD(3);
t37 = t50 - t75;
t60 = t31 * t57 + t32 * t55 + (-t37 * t57 - t38 * t55) * qJD(3);
t36 = t57 * t70;
t16 = t18 * t51;
t13 = t56 * t26 - t92;
t12 = -t54 * t26 - t90;
t10 = t56 * t25 - t92;
t8 = -t28 ^ 2 + t30 ^ 2;
t5 = -t85 + (t28 - t74) * t51;
t4 = -t20 * qJD(4) - t54 * t35 + t56 * t36;
t3 = t19 * qJD(4) + t56 * t35 + t54 * t36;
t1 = t67 + t98;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t86, 0, t31 * t55 - t32 * t57 + (-t37 * t55 + t38 * t57) * qJD(3), 0, 0, 0, 0, 0, 0, -t16, t95, t62 + t96, t1 * t34 - t10 * t18 - t11 * t17 - t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t64, t84 * t99, t86, -0.2e1 * t64, -t87, 0, -pkin(5) * t86 + t55 * t65, pkin(5) * t87 + t57 * t65, t60, t60 * pkin(5), -t14 * t34 - t30 * t17, -t62 + t96, -t95, t15 * t33 + t28 * t18, -t16, 0, t49 * t15 + t39 * t18 + t4 * t51 + (qJD(2) * t33 + t28) * t76, t30 * t100 - t49 * t14 - t39 * t17 - t3 * t51, -t1 * t33 + t10 * t17 - t11 * t18 + t19 * t14 - t20 * t15 - t2 * t34 - t3 * t28 - t4 * t30, t1 * t20 + t10 * t4 + t39 * t100 + t11 * t3 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t84 * t59, 0, t78, 0, 0, t59 * pkin(2) * t55, pkin(2) * t88 - t47 + (t37 + t75) * qJD(3), 0, 0, t94, t8, t5, -t94, 0, 0, -t28 * t77 - t12 * t51 - t93 + (t72 * t54 - t90) * qJD(4) + t68, -t30 * t77 + t13 * t51 + (t72 * qJD(4) - t23) * t56 + t61, (t11 + t12) * t30 + (-t10 + t13) * t28 + (t14 * t56 - t15 * t54 + (-t28 * t56 + t30 * t54) * qJD(4)) * pkin(3), -t10 * t12 - t11 * t13 + (-t39 * t82 + t1 * t54 + t2 * t56 + (-t10 * t54 + t11 * t56) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t8, t5, -t94, 0, 0, t11 * t51 + t2 - t93, t10 * t51 + t61 - t98, 0, 0;];
tauc_reg = t6;
