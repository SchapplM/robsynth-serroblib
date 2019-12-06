% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:22
% EndTime: 2019-12-05 15:03:24
% DurationCPUTime: 0.29s
% Computational Cost: add. (382->67), mult. (921->97), div. (0->0), fcn. (658->6), ass. (0->65)
t58 = cos(pkin(8));
t68 = cos(qJ(3));
t37 = t68 * t58;
t25 = sin(pkin(8));
t27 = sin(qJ(3));
t65 = t27 * t25;
t15 = -t37 + t65;
t46 = t27 * t58;
t51 = qJD(1) * qJD(3);
t47 = t25 * t51;
t10 = t46 * t51 + t68 * t47;
t16 = t68 * t25 + t46;
t12 = t16 * qJD(1);
t52 = qJD(3) * qJ(4);
t9 = t12 + t52;
t59 = t9 * qJD(3);
t70 = -t59 + t10;
t20 = qJD(1) * t37;
t48 = qJD(1) * t65;
t11 = -t20 + t48;
t19 = qJD(3) * t20;
t69 = (-t11 + t48) * qJD(3) - t19;
t26 = sin(qJ(5));
t28 = cos(qJ(5));
t29 = -pkin(3) - pkin(6);
t35 = qJD(4) + t11;
t5 = t29 * qJD(3) + t35;
t3 = -t26 * qJD(2) + t28 * t5;
t1 = t3 * qJD(5) + t26 * t10;
t4 = t28 * qJD(2) + t26 * t5;
t7 = t28 * t10;
t2 = -t4 * qJD(5) + t7;
t32 = -(t26 * t3 - t28 * t4) * qJD(5) + t1 * t26 + t2 * t28;
t67 = t10 * t15;
t66 = t26 * t28;
t30 = qJD(5) ^ 2;
t64 = t30 * t26;
t63 = t30 * t28;
t23 = t26 ^ 2;
t24 = t28 ^ 2;
t62 = t23 - t24;
t61 = t23 + t24;
t31 = qJD(3) ^ 2;
t60 = -t30 - t31;
t57 = qJD(5) * t26;
t56 = qJD(5) * t28;
t55 = t12 * qJD(3);
t13 = t15 * qJD(3);
t54 = t13 * qJD(3);
t14 = t16 * qJD(3);
t53 = t14 * qJD(3);
t50 = qJD(3) * qJD(5);
t49 = t31 * t66;
t45 = t50 * t66;
t22 = qJD(3) * qJD(4);
t36 = -t27 * t47 + t19;
t6 = -t22 - t36;
t42 = -t9 * t13 - t6 * t16;
t41 = t26 * t4 + t28 * t3;
t39 = -t12 + t9 + t52;
t38 = t55 - t10;
t34 = -t6 * qJ(4) + t35 * t9;
t33 = t35 * qJD(3) - t29 * t30 - t6;
t8 = -qJD(3) * pkin(3) + t35;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t54, 0, t11 * t14 - t12 * t13 + t36 * t16 + t67, 0, 0, 0, 0, 0, 0, 0, t53, -t54, t8 * t14 + t42 + t67, 0, 0, 0, 0, 0, 0, t14 * t56 - t15 * t64 + (-t13 * t26 + t16 * t56) * qJD(3), -t14 * t57 - t15 * t63 + (-t13 * t28 - t16 * t57) * qJD(3), -t61 * t53, t41 * t14 + t32 * t15 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t64, 0, -t41 * qJD(5) + t1 * t28 - t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0.2e1 * t22 - t69, -t10 * pkin(3) - t8 * t12 + t34, -0.2e1 * t45, 0.2e1 * t62 * t50, -t64, 0.2e1 * t45, -t63, 0, t33 * t26 + t39 * t56, t33 * t28 - t39 * t57, t61 * t55 - t32, -t41 * t12 + t32 * t29 + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t70, 0, 0, 0, 0, 0, 0, t60 * t26, t60 * t28, 0, t32 - t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t62 * t31, 0, -t49, 0, 0, -t28 * t59 + t7, -t70 * t26, 0, 0;];
tauc_reg = t17;
