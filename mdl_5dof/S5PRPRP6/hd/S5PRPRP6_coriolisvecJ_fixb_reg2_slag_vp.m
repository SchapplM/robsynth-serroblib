% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:19
% EndTime: 2019-12-05 15:41:21
% DurationCPUTime: 0.41s
% Computational Cost: add. (357->95), mult. (787->123), div. (0->0), fcn. (375->4), ass. (0->76)
t31 = sin(qJ(4));
t32 = sin(qJ(2));
t61 = t32 * qJD(1);
t53 = qJD(2) * t61;
t25 = t31 * t53;
t34 = cos(qJ(2));
t60 = t34 * qJD(1);
t49 = qJD(3) - t60;
t35 = -pkin(2) - pkin(6);
t80 = qJD(2) * t35;
t13 = t49 + t80;
t33 = cos(qJ(4));
t76 = t33 * t13;
t2 = t25 + (qJD(5) + t76) * qJD(4);
t27 = t33 * t53;
t63 = qJD(4) * t31;
t3 = t13 * t63 - t27;
t50 = qJD(4) * pkin(4) - qJD(5);
t6 = -t50 - t76;
t58 = qJD(4) * qJ(5);
t8 = t31 * t13 + t58;
t45 = t31 * t6 + t33 * t8;
t39 = t45 * qJD(4) + t2 * t31 - t3 * t33;
t15 = t31 * pkin(4) - t33 * qJ(5) + qJ(3);
t64 = qJD(2) * t15;
t7 = t61 + t64;
t65 = t7 * qJD(2);
t82 = t39 - t65;
t29 = t31 ^ 2;
t30 = t33 ^ 2;
t68 = t29 + t30;
t81 = t32 * t68;
t44 = pkin(4) * t33 + qJ(5) * t31;
t9 = t44 * qJD(4) - t33 * qJD(5) + qJD(3);
t1 = (t9 + t60) * qJD(2);
t36 = qJD(4) ^ 2;
t75 = t35 * t36;
t79 = (-t9 + t60) * qJD(2) - t1 + t75;
t16 = (qJD(3) + t60) * qJD(2);
t78 = t16 * t32;
t59 = qJD(2) * qJ(3);
t23 = t59 + t61;
t77 = t23 * t34;
t74 = t36 * t31;
t73 = t36 * t33;
t37 = qJD(2) ^ 2;
t72 = t37 * t32;
t71 = t37 * t34;
t70 = t68 * t53;
t69 = t29 - t30;
t67 = t36 + t37;
t66 = qJD(2) * pkin(2);
t62 = t23 * qJD(2);
t57 = qJD(2) * qJD(4);
t55 = t67 * t34;
t54 = t33 * t57;
t51 = -t23 + t61;
t48 = t31 * t54;
t46 = t31 * t8 - t33 * t6;
t43 = t16 * qJ(3) + t23 * qJD(3);
t42 = -t51 + t59;
t41 = qJD(4) * (-t61 + t7 + t64);
t40 = t49 * qJD(2) + t16 - t75;
t28 = t33 * t37 * t31;
t22 = -0.2e1 * t48;
t21 = 0.2e1 * t48;
t20 = t49 - t66;
t19 = t67 * t33;
t18 = t67 * t31;
t17 = t69 * t37;
t14 = t44 * qJD(2);
t12 = t69 * t57;
t11 = t68 * t72;
t5 = -0.2e1 * t31 * t32 * t57 + t33 * t55;
t4 = t31 * t55 + 0.2e1 * t32 * t54;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t71, t78 + (t77 + (t20 - t60) * t32) * qJD(2), 0, 0, 0, 0, 0, 0, t4, t5, -t11, t78 + (t77 + (t13 - t60) * t81) * qJD(2), 0, 0, 0, 0, 0, 0, t4, -t11, -t5, (t46 * qJD(2) + t1) * t32 - t82 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), (-t77 + (-t20 - t66) * t32) * qJD(1) + t43, t22, 0.2e1 * t12, -t74, t21, -t73, 0, t42 * t33 * qJD(4) + t40 * t31, t40 * t33 - t42 * t63, 0, (-t77 + (-t13 + t80) * t81) * qJD(1) + t43, t22, -t74, -0.2e1 * t12, 0, t73, t21, -t79 * t31 + t33 * t41, -t39 + t70, t31 * t41 + t79 * t33, t1 * t15 + t7 * t9 + (-t46 * t32 - t34 * t7) * qJD(1) + t39 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t51 * qJD(2), 0, 0, 0, 0, 0, 0, -t18, -t19, 0, -t62 + t70, 0, 0, 0, 0, 0, 0, -t18, 0, t19, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t17, 0, -t28, 0, 0, -t33 * t62 + t27, t31 * t62 - t25, 0, 0, t28, 0, t17, 0, 0, -t28, t27 + (-t14 * t31 - t33 * t7) * qJD(2), ((t8 - t58) * t33 + (t50 + t6) * t31) * qJD(2), 0.2e1 * qJD(4) * qJD(5) + t25 + (t14 * t33 - t31 * t7) * qJD(2), -t3 * pkin(4) + t2 * qJ(5) + t8 * qJD(5) - t45 * t13 - t7 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t30 * t37 - t36, -t8 * qJD(4) + t33 * t65 + t3;];
tauc_reg = t10;
