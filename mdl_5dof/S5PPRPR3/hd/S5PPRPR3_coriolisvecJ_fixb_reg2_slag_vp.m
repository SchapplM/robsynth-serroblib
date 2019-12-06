% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:22
% EndTime: 2019-12-05 15:05:24
% DurationCPUTime: 0.40s
% Computational Cost: add. (604->89), mult. (1553->159), div. (0->0), fcn. (1243->8), ass. (0->71)
t46 = sin(qJ(3));
t48 = cos(qJ(3));
t42 = sin(pkin(8));
t67 = qJD(1) * t42;
t33 = t46 * qJD(2) + t48 * t67;
t43 = cos(pkin(9));
t23 = t43 * t33;
t41 = sin(pkin(9));
t80 = t48 * qJD(2) - t46 * t67;
t13 = t41 * t80 + t23;
t28 = t33 * qJD(3);
t29 = t80 * qJD(3);
t9 = t43 * t28 + t41 * t29;
t81 = t13 * qJD(3) - t9;
t30 = t41 * t46 - t43 * t48;
t31 = t41 * t48 + t43 * t46;
t21 = t31 * t42;
t79 = t9 * t21;
t78 = t9 * t30;
t77 = t41 * t33;
t45 = sin(qJ(5));
t47 = cos(qJ(5));
t74 = t45 * t47;
t49 = qJD(5) ^ 2;
t73 = t49 * t45;
t72 = t49 * t47;
t50 = qJD(3) ^ 2;
t71 = t50 * t46;
t70 = t50 * t48;
t25 = qJD(3) * pkin(3) + t80;
t12 = t41 * t25 + t23;
t39 = t45 ^ 2;
t40 = t47 ^ 2;
t69 = t39 - t40;
t68 = t39 + t40;
t66 = qJD(3) * t42;
t65 = qJD(5) * t45;
t64 = qJD(5) * t47;
t14 = t43 * t80 - t77;
t62 = t14 * qJD(3);
t27 = t30 * qJD(3);
t61 = t27 * qJD(3);
t60 = qJD(3) * qJD(5);
t59 = t50 * t74;
t10 = -t41 * t28 + t43 * t29;
t11 = t43 * t25 - t77;
t7 = -qJD(3) * pkin(4) - t11;
t57 = -qJD(3) * t7 - t10;
t56 = t60 * t74;
t44 = cos(pkin(8));
t35 = -t44 * qJD(1) + qJD(4);
t8 = qJD(3) * pkin(6) + t12;
t5 = t47 * t35 - t45 * t8;
t6 = t45 * t35 + t47 * t8;
t55 = t45 * t5 - t47 * t6;
t22 = t30 * t42;
t16 = -t47 * t22 - t44 * t45;
t15 = t45 * t22 - t44 * t47;
t36 = t41 * pkin(3) + pkin(6);
t54 = -t36 * t49 + t81;
t37 = -t43 * pkin(3) - pkin(4);
t53 = qJD(5) * (qJD(3) * t37 + t14 + t7);
t1 = t5 * qJD(5) + t47 * t10;
t2 = -t6 * qJD(5) - t45 * t10;
t51 = t1 * t47 - t2 * t45 + (-t45 * t6 - t47 * t5) * qJD(5);
t26 = t31 * qJD(3);
t18 = t31 * t66;
t17 = t30 * t66;
t4 = -qJD(5) * t16 + t45 * t18;
t3 = qJD(5) * t15 - t47 * t18;
t19 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42 * t70, t42 * t71, 0, (t28 * t46 + t29 * t48 + (-t33 * t46 - t48 * t80) * qJD(3)) * t42, 0, 0, 0, 0, 0, 0, t17 * qJD(3), t18 * qJD(3), 0, -t10 * t22 + t11 * t17 - t12 * t18 + t79, 0, 0, 0, 0, 0, 0, t4 * qJD(5) + (t17 * t47 + t21 * t65) * qJD(3), -t3 * qJD(5) + (-t17 * t45 + t21 * t64) * qJD(3), (t3 * t47 - t4 * t45 + (-t15 * t47 - t16 * t45) * qJD(5)) * qJD(3), t1 * t16 + t2 * t15 - t7 * t17 + t6 * t3 + t5 * t4 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t70, 0, -t28 * t48 + t29 * t46 + (t33 * t48 - t46 * t80) * qJD(3), 0, 0, 0, 0, 0, 0, -t26 * qJD(3), t61, 0, t10 * t31 - t11 * t26 - t12 * t27 + t78, 0, 0, 0, 0, 0, 0, t27 * t65 - t31 * t72 + (-t26 * t47 + t30 * t65) * qJD(3), t27 * t64 + t31 * t73 + (t26 * t45 + t30 * t64) * qJD(3), -t68 * t61, t7 * t26 + t55 * t27 + t51 * t31 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t10 + t62, 0, t11 * t13 - t12 * t14 + (t10 * t41 - t43 * t9) * pkin(3), 0.2e1 * t56, -0.2e1 * t69 * t60, t72, -0.2e1 * t56, -t73, 0, t45 * t53 + t54 * t47, -t54 * t45 + t47 * t53, -t68 * t62 + t51, -t7 * t13 + t55 * t14 + t51 * t36 + t9 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t72, 0, -t55 * qJD(5) + t1 * t45 + t2 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t69 * t50, 0, t59, 0, 0, t57 * t45, t57 * t47, 0, 0;];
tauc_reg = t19;
