% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:12
% EndTime: 2020-01-03 11:59:13
% DurationCPUTime: 0.35s
% Computational Cost: add. (672->105), mult. (1330->153), div. (0->0), fcn. (735->6), ass. (0->85)
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t50 = qJD(1) + qJD(2);
t58 = cos(qJ(2));
t87 = pkin(1) * qJD(1);
t78 = t58 * t87;
t35 = t50 * pkin(2) + t78;
t54 = cos(pkin(8));
t56 = sin(qJ(2));
t80 = t56 * t87;
t39 = t54 * t80;
t53 = sin(pkin(8));
t18 = t53 * t35 + t39;
t72 = t18 + (pkin(7) + qJ(5)) * t50;
t7 = t57 * qJD(3) - t72 * t55;
t8 = t55 * qJD(3) + t72 * t57;
t85 = qJD(4) * pkin(4);
t6 = t7 + t85;
t97 = t6 - t7;
t96 = t54 * pkin(2);
t95 = t57 * pkin(4);
t94 = t53 * t56;
t93 = t54 * t56;
t59 = qJD(4) ^ 2;
t92 = t59 * t55;
t38 = t53 * t80;
t17 = t54 * t35 - t38;
t13 = -t50 * pkin(3) - t17;
t86 = pkin(1) * qJD(2);
t28 = (t53 * t58 + t93) * t86;
t22 = qJD(1) * t28;
t91 = t13 * qJD(4) * t57 + t22 * t55;
t44 = t58 * pkin(1) + pkin(2);
t90 = pkin(1) * t93 + t53 * t44;
t51 = t55 ^ 2;
t52 = t57 ^ 2;
t89 = -t51 - t52;
t88 = t51 - t52;
t26 = pkin(7) + t90;
t84 = -qJ(5) - t26;
t42 = t53 * pkin(2) + pkin(7);
t83 = -qJ(5) - t42;
t81 = t55 * qJD(4);
t79 = pkin(4) * t81;
t77 = t50 * t81;
t76 = -pkin(3) - t95;
t30 = (t54 * t58 - t94) * t86;
t23 = qJD(1) * t30;
t71 = qJD(5) * t50 + t23;
t2 = t7 * qJD(4) + t71 * t57;
t3 = -qJD(4) * t8 - t71 * t55;
t75 = t2 * t57 - t3 * t55;
t74 = -t13 * t50 - t23;
t73 = -pkin(1) * t94 + t54 * t44;
t70 = qJD(4) * t84;
t69 = qJD(4) * t83;
t25 = -pkin(3) - t73;
t68 = -t55 * t8 - t57 * t6;
t67 = t6 * t55 - t8 * t57;
t66 = (-qJD(2) + t50) * t87;
t65 = (-qJD(1) - t50) * t86;
t64 = t26 * t59 + t28 * t50;
t27 = t53 * t78 + t39;
t63 = -t27 * t50 + t42 * t59;
t61 = qJD(4) * (t25 * t50 - t30);
t29 = t54 * t78 - t38;
t60 = qJD(4) * ((-pkin(3) - t96) * t50 + t29);
t12 = pkin(4) * t77 + t22;
t49 = t50 ^ 2;
t48 = t57 * qJ(5);
t47 = t59 * t57;
t45 = t57 * qJD(5);
t34 = 0.2e1 * t57 * t77;
t33 = t57 * t42 + t48;
t32 = t83 * t55;
t24 = -0.2e1 * t88 * t50 * qJD(4);
t21 = -t55 * qJD(5) + t57 * t69;
t20 = t55 * t69 + t45;
t16 = t57 * t26 + t48;
t15 = t84 * t55;
t10 = t13 * t81;
t9 = t76 * t50 + qJD(5) - t17;
t5 = (-qJD(5) - t30) * t55 + t57 * t70;
t4 = t57 * t30 + t55 * t70 + t45;
t1 = [0, 0, 0, 0, t56 * t65, t58 * t65, -t17 * t28 + t18 * t30 - t22 * t73 + t23 * t90, t34, t24, t47, -t92, 0, t10 + t55 * t61 + (-t22 - t64) * t57, t64 * t55 + t57 * t61 + t91, (t4 * t57 - t5 * t55) * t50 + ((-t15 * t57 - t16 * t55) * t50 + t68) * qJD(4) + t75, t2 * t16 + t8 * t4 + t3 * t15 + t6 * t5 + t12 * (t25 - t95) + t9 * (t28 + t79); 0, 0, 0, 0, t56 * t66, t58 * t66, t17 * t27 - t18 * t29 + (-t22 * t54 + t23 * t53) * pkin(2), t34, t24, t47, -t92, 0, t10 + t55 * t60 + (-t22 - t63) * t57, t63 * t55 + t57 * t60 + t91, t68 * qJD(4) + (t20 * t57 - t21 * t55 + t89 * t29 + (-t32 * t57 - t33 * t55) * qJD(4)) * t50 + t75, t2 * t33 + t8 * t20 + t3 * t32 + t6 * t21 + t12 * (t76 - t96) + (-t27 + t79) * t9 + t67 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t47, 0, -t67 * qJD(4) + t2 * t55 + t3 * t57; 0, 0, 0, 0, 0, 0, 0, -t55 * t49 * t57, t88 * t49, 0, 0, 0, t74 * t55, t74 * t57, (-t85 + t97) * t57 * t50, t97 * t8 + (-t50 * t55 * t9 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t49, t67 * t50 + t12;];
tauc_reg = t1;
