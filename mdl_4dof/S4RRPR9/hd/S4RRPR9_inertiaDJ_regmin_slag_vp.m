% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:04
% EndTime: 2019-12-31 17:10:06
% DurationCPUTime: 0.41s
% Computational Cost: add. (329->81), mult. (955->184), div. (0->0), fcn. (795->6), ass. (0->59)
t39 = sin(pkin(7));
t40 = cos(pkin(7));
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t68 = -t41 * t39 + t43 * t40;
t35 = -t40 * pkin(3) - pkin(2);
t67 = 0.2e1 * t35;
t44 = cos(qJ(2));
t66 = t39 * t44;
t42 = sin(qJ(2));
t65 = t40 * t42;
t64 = t40 * t44;
t61 = pkin(6) + qJ(3);
t19 = -t42 * qJD(3) + (pkin(2) * t42 - qJ(3) * t44) * qJD(2);
t58 = t42 * qJD(2);
t54 = pkin(5) * t58;
t10 = t40 * t19 + t39 * t54;
t48 = -t44 * pkin(2) - t42 * qJ(3);
t28 = -pkin(1) + t48;
t33 = pkin(5) * t64;
t14 = t39 * t28 + t33;
t60 = qJD(3) * t44;
t59 = qJD(4) * t42;
t57 = t44 * qJD(2);
t56 = pkin(5) * t66;
t55 = -0.2e1 * pkin(1) * qJD(2);
t36 = pkin(5) * t57;
t53 = t39 * t57;
t52 = t42 * t57;
t51 = 0.2e1 * (t39 ^ 2 + t40 ^ 2) * qJD(3);
t12 = -t39 * t42 * pkin(6) + t14;
t24 = t40 * t28;
t9 = -pkin(6) * t65 + t24 + (-pkin(5) * t39 - pkin(3)) * t44;
t50 = t43 * t12 + t41 * t9;
t49 = t41 * t12 - t43 * t9;
t15 = t39 * t19;
t11 = -t40 * t54 + t15;
t47 = -t10 * t39 + t11 * t40;
t30 = t61 * t39;
t31 = t61 * t40;
t46 = -t43 * t30 - t41 * t31;
t45 = -t41 * t30 + t43 * t31;
t26 = t43 * t39 + t41 * t40;
t27 = (pkin(3) * t39 + pkin(5)) * t42;
t22 = pkin(3) * t53 + t36;
t21 = t26 * qJD(4);
t20 = t68 * qJD(4);
t18 = t68 * t42;
t17 = t26 * t42;
t13 = t24 - t56;
t8 = t15 + (-pkin(5) * t65 - pkin(6) * t66) * qJD(2);
t7 = t26 * t57 + t68 * t59;
t6 = -t26 * t59 + t57 * t68;
t5 = (pkin(3) * t42 - pkin(6) * t64) * qJD(2) + t10;
t4 = -t26 * qJD(3) - t45 * qJD(4);
t3 = -qJD(3) * t68 - t46 * qJD(4);
t2 = -t50 * qJD(4) - t41 * t8 + t43 * t5;
t1 = t49 * qJD(4) - t41 * t5 - t43 * t8;
t16 = [0, 0, 0, 0.2e1 * t52, 0.2e1 * (-t42 ^ 2 + t44 ^ 2) * qJD(2), 0, 0, 0, t42 * t55, t44 * t55, -0.2e1 * t10 * t44 + 0.2e1 * (t13 + 0.2e1 * t56) * t58, 0.2e1 * t11 * t44 + 0.2e1 * (-t14 + 0.2e1 * t33) * t58, 0.2e1 * (-t10 * t40 - t11 * t39) * t42 + 0.2e1 * (-t13 * t40 - t14 * t39) * t57, 0.2e1 * pkin(5) ^ 2 * t52 + 0.2e1 * t13 * t10 + 0.2e1 * t14 * t11, 0.2e1 * t18 * t6, -0.2e1 * t6 * t17 - 0.2e1 * t18 * t7, 0.2e1 * t18 * t58 - 0.2e1 * t6 * t44, -0.2e1 * t17 * t58 + 0.2e1 * t7 * t44, -0.2e1 * t52, 0.2e1 * t22 * t17 - 0.2e1 * t2 * t44 + 0.2e1 * t27 * t7 - 0.2e1 * t49 * t58, -0.2e1 * t1 * t44 + 0.2e1 * t22 * t18 + 0.2e1 * t27 * t6 - 0.2e1 * t50 * t58; 0, 0, 0, 0, 0, t57, -t58, 0, -t36, t54, t39 * t60 + (t48 * t39 - t33) * qJD(2), t40 * t60 + (t48 * t40 + t56) * qJD(2), t47, -pkin(2) * t36 + (-t13 * t39 + t14 * t40) * qJD(3) + t47 * qJ(3), t18 * t20 + t6 * t26, -t20 * t17 - t18 * t21 - t26 * t7 + t6 * t68, -t20 * t44 + t26 * t58, t21 * t44 + t58 * t68, 0, t27 * t21 - t22 * t68 + t35 * t7 - t4 * t44 + t46 * t58, t27 * t20 + t22 * t26 - t3 * t44 + t35 * t6 - t45 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, qJ(3) * t51, 0.2e1 * t26 * t20, 0.2e1 * t20 * t68 - 0.2e1 * t26 * t21, 0, 0, 0, t21 * t67, t20 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t40 * t57, 0, t36, 0, 0, 0, 0, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t7, t58, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t21, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t16;
