% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:13
% EndTime: 2019-12-05 18:12:15
% DurationCPUTime: 0.45s
% Computational Cost: add. (1255->75), mult. (2828->139), div. (0->0), fcn. (2920->8), ass. (0->62)
t42 = sin(pkin(9));
t43 = cos(pkin(9));
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t33 = t49 * t42 + t46 * t43;
t31 = t33 * qJD(3);
t52 = t46 * t42 - t49 * t43;
t69 = pkin(6) + qJ(2);
t34 = t69 * t42;
t35 = t69 * t43;
t54 = -t49 * t34 - t46 * t35;
t72 = t52 * qJD(2) - t54 * qJD(3);
t18 = -t31 * pkin(7) - t72;
t30 = t52 * qJD(3);
t53 = t46 * t34 - t49 * t35;
t50 = -t33 * qJD(2) + t53 * qJD(3);
t19 = t30 * pkin(7) + t50;
t45 = sin(qJ(4));
t48 = cos(qJ(4));
t21 = -t33 * pkin(7) + t54;
t22 = -pkin(7) * t52 - t53;
t57 = t48 * t21 - t45 * t22;
t7 = -t57 * qJD(4) - t48 * t18 - t45 * t19;
t38 = -t43 * pkin(2) - pkin(1);
t71 = 0.2e1 * t38;
t70 = t31 * pkin(3);
t44 = sin(qJ(5));
t67 = pkin(3) * qJD(4);
t64 = t45 * t67;
t66 = qJD(5) * t44;
t68 = t45 * pkin(3) * t66 + t44 * t64;
t47 = cos(qJ(5));
t65 = qJD(5) * t47;
t63 = t48 * t67;
t62 = pkin(4) * t66;
t61 = pkin(4) * t65;
t39 = t48 * pkin(3) + pkin(4);
t60 = (-pkin(4) - t39) * qJD(5);
t59 = 0.2e1 * (t42 ^ 2 + t43 ^ 2) * qJD(2);
t56 = -t45 * t21 - t48 * t22;
t26 = t48 * t33 - t45 * t52;
t55 = -t45 * t33 - t48 * t52;
t15 = t44 * t26 - t47 * t55;
t16 = t47 * t26 + t44 * t55;
t27 = pkin(3) * t52 + t38;
t8 = t56 * qJD(4) - t45 * t18 + t48 * t19;
t51 = (-t45 * t65 + (-t44 * t48 - t45 * t47) * qJD(4)) * pkin(3);
t24 = -t39 * t66 + t51;
t23 = (-qJD(5) * t39 - t63) * t47 + t68;
t20 = -pkin(4) * t55 + t27;
t14 = t26 * qJD(4) - t45 * t30 + t48 * t31;
t13 = t55 * qJD(4) - t48 * t30 - t45 * t31;
t11 = t14 * pkin(4) + t70;
t10 = pkin(8) * t55 - t56;
t9 = -t26 * pkin(8) + t57;
t6 = t16 * qJD(5) + t44 * t13 + t47 * t14;
t5 = -t15 * qJD(5) + t47 * t13 - t44 * t14;
t4 = -t13 * pkin(8) + t8;
t3 = -t14 * pkin(8) - t7;
t2 = -t44 * t3 + t47 * t4 + (-t10 * t47 - t44 * t9) * qJD(5);
t1 = -t47 * t3 - t44 * t4 + (t10 * t44 - t47 * t9) * qJD(5);
t12 = [0, 0, 0, 0, 0, t59, qJ(2) * t59, -0.2e1 * t33 * t30, 0.2e1 * t30 * t52 - 0.2e1 * t33 * t31, 0, 0, 0, t31 * t71, -t30 * t71, 0.2e1 * t26 * t13, 0.2e1 * t13 * t55 - 0.2e1 * t26 * t14, 0, 0, 0, 0.2e1 * t27 * t14 - 0.2e1 * t55 * t70, 0.2e1 * t27 * t13 + 0.2e1 * t26 * t70, 0.2e1 * t16 * t5, -0.2e1 * t5 * t15 - 0.2e1 * t16 * t6, 0, 0, 0, 0.2e1 * t11 * t15 + 0.2e1 * t20 * t6, 0.2e1 * t11 * t16 + 0.2e1 * t20 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, 0, 0, 0, 0, t14, t13, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t31, 0, t50, t72, 0, 0, t13, -t14, 0, t8, t7, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t64, -0.2e1 * t63, 0, 0, 0, 0, 0, 0.2e1 * t24, 0.2e1 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, 0, t8, t7, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63, 0, 0, 0, 0, 0, t44 * t60 + t51, (t60 - t63) * t47 + t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t62, -0.2e1 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
