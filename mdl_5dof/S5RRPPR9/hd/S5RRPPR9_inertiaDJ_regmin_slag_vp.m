% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:47
% EndTime: 2019-12-31 19:41:49
% DurationCPUTime: 0.47s
% Computational Cost: add. (238->79), mult. (574->160), div. (0->0), fcn. (370->4), ass. (0->59)
t67 = pkin(6) - qJ(4);
t34 = sin(qJ(2));
t22 = t34 * qJ(3);
t36 = cos(qJ(2));
t66 = -t36 * pkin(2) - t22;
t35 = cos(qJ(5));
t30 = t35 ^ 2;
t33 = sin(qJ(5));
t62 = t33 ^ 2 - t30;
t46 = t62 * qJD(5);
t32 = qJ(3) + pkin(4);
t56 = t34 * qJD(2);
t52 = pkin(6) * t56;
t54 = qJ(4) * qJD(2);
t6 = t36 * qJD(4) - t34 * t54 + t52;
t37 = -pkin(2) - pkin(3);
t27 = -pkin(7) + t37;
t64 = t27 * t34;
t65 = (t32 * t36 + t64) * qJD(5) + t6;
t38 = 2 * qJD(3);
t20 = t36 * qJD(2);
t63 = qJ(3) * t20 + t34 * qJD(3);
t31 = t36 ^ 2;
t61 = t34 ^ 2 - t31;
t60 = qJD(2) * t35;
t21 = qJD(5) * t33;
t59 = qJD(5) * t35;
t58 = qJD(5) * t36;
t14 = t67 * t36;
t57 = t14 * qJD(5);
t55 = t36 * qJD(3);
t53 = -0.2e1 * pkin(1) * qJD(2);
t12 = -pkin(1) + t66;
t51 = t33 * t58;
t50 = t35 * t58;
t49 = t33 * t59;
t48 = t34 * t20;
t47 = t35 * t56;
t11 = t36 * pkin(3) - t12;
t45 = t61 * qJD(2);
t44 = t33 * t47;
t13 = t67 * t34;
t4 = t34 * pkin(4) + t36 * pkin(7) + t11;
t43 = t35 * t13 + t33 * t4;
t42 = t33 * t13 - t35 * t4;
t40 = t66 * qJD(2) + t55;
t39 = -t55 - t57 + (-t27 * t36 + t32 * t34) * qJD(2);
t26 = qJ(3) * t38;
t18 = pkin(6) * t20;
t15 = 0.2e1 * t48;
t10 = -t33 * t20 - t34 * t59;
t9 = -t35 * t20 + t34 * t21;
t8 = pkin(2) * t56 - t63;
t7 = -t34 * qJD(4) - t36 * t54 + t18;
t5 = t37 * t56 + t63;
t3 = (pkin(4) * t36 + t64) * qJD(2) + t63;
t2 = -t43 * qJD(5) + t35 * t3 - t33 * t7;
t1 = t42 * qJD(5) - t33 * t3 - t35 * t7;
t16 = [0, 0, 0, t15, -0.2e1 * t45, 0, 0, 0, t34 * t53, t36 * t53, 0.2e1 * t12 * t56 - 0.2e1 * t8 * t36, 0, -0.2e1 * t12 * t20 - 0.2e1 * t8 * t34, 0.2e1 * t12 * t8, 0.2e1 * t11 * t20 + 0.2e1 * t5 * t34, 0.2e1 * t11 * t56 - 0.2e1 * t5 * t36, -0.2e1 * t7 * t34 + 0.2e1 * t6 * t36 + 0.2e1 * (-t13 * t36 + t14 * t34) * qJD(2), 0.2e1 * t11 * t5 + 0.2e1 * t13 * t7 - 0.2e1 * t14 * t6, -0.2e1 * t30 * t48 - 0.2e1 * t31 * t49, 0.2e1 * t31 * t46 + 0.4e1 * t36 * t44, 0.2e1 * t34 * t51 + 0.2e1 * t61 * t60, -0.2e1 * t33 * t45 + 0.2e1 * t34 * t50, t15, 0.2e1 * (t14 * t33 * qJD(2) + t2) * t34 + 0.2e1 * (-t42 * qJD(2) + t6 * t33 - t35 * t57) * t36, 0.2e1 * (t14 * t60 + t1) * t34 + 0.2e1 * (-t43 * qJD(2) + t33 * t57 + t6 * t35) * t36; 0, 0, 0, 0, 0, t20, -t56, 0, -t18, t52, -t18, t40, -t52, t40 * pkin(6), -t6, t7, -t55 + (-t36 * t37 + t22) * qJD(2), -t6 * qJ(3) + t14 * qJD(3) + t7 * t37, -t36 * t46 - t44, -0.4e1 * t36 * t49 + t62 * t56, t10, t9, 0, t39 * t33 - t65 * t35, t65 * t33 + t39 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t26, t38, 0, 0, t26, 0.2e1 * t49, -0.2e1 * t46, 0, 0, 0, 0.2e1 * qJD(3) * t35 - 0.2e1 * t32 * t21, -0.2e1 * qJD(3) * t33 - 0.2e1 * t32 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, t18, 0, 0, -t20, t7, 0, 0, 0, 0, 0, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t56, 0, t5, 0, 0, 0, 0, 0, -t9, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 + t51, -t33 * t56 + t50, t20, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t21, 0, -t27 * t59, t27 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t16;
