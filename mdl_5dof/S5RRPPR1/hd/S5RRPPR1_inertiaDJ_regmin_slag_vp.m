% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:27
% EndTime: 2019-12-05 18:18:29
% DurationCPUTime: 0.29s
% Computational Cost: add. (189->48), mult. (536->86), div. (0->0), fcn. (415->8), ass. (0->47)
t35 = sin(pkin(9));
t37 = cos(pkin(9));
t51 = t35 ^ 2 + t37 ^ 2;
t66 = qJD(4) * t51;
t67 = 0.2e1 * t66;
t36 = sin(pkin(8));
t38 = cos(pkin(8));
t42 = cos(qJ(2));
t50 = pkin(1) * qJD(2);
t48 = t42 * t50;
t40 = sin(qJ(2));
t49 = t40 * t50;
t15 = -t36 * t49 + t38 * t48;
t12 = qJD(4) + t15;
t65 = t12 * t51;
t39 = sin(qJ(5));
t41 = cos(qJ(5));
t20 = t39 * t35 - t41 * t37;
t29 = t42 * pkin(1) + pkin(2);
t45 = -t36 * t40 * pkin(1) + t38 * t29;
t43 = -pkin(3) - t45;
t62 = t37 * pkin(4);
t11 = t43 - t62;
t56 = t38 * t40;
t14 = (t36 * t42 + t56) * t50;
t21 = t41 * t35 + t39 * t37;
t17 = t21 * qJD(5);
t64 = t11 * t17 + t14 * t20;
t16 = t20 * qJD(5);
t63 = -t11 * t16 + t14 * t21;
t60 = t14 * t35;
t59 = t14 * t37;
t47 = -t38 * pkin(2) - pkin(3);
t22 = t47 - t62;
t58 = t22 * t16;
t57 = t22 * t17;
t53 = pkin(1) * t56 + t36 * t29;
t32 = t37 * pkin(7);
t28 = t36 * pkin(2) + qJ(4);
t19 = t37 * t28 + t32;
t18 = (-pkin(7) - t28) * t35;
t13 = qJ(4) + t53;
t8 = t37 * t13 + t32;
t7 = (-pkin(7) - t13) * t35;
t6 = -0.2e1 * t21 * t16;
t1 = 0.2e1 * t16 * t20 - 0.2e1 * t21 * t17;
t2 = [0, 0, 0, 0, -0.2e1 * t49, -0.2e1 * t48, -0.2e1 * t45 * t14 + 0.2e1 * t53 * t15, -0.2e1 * t59, 0.2e1 * t60, 0.2e1 * t65, 0.2e1 * t13 * t65 + 0.2e1 * t43 * t14, t6, t1, 0, 0, 0, 0.2e1 * t64, 0.2e1 * t63; 0, 0, 0, 0, -t49, -t48, (-t14 * t38 + t15 * t36) * pkin(2), -t59, t60, t66 + t65, t13 * t66 + t14 * t47 + t28 * t65, t6, t1, 0, 0, 0, t57 + t64, -t58 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t28 * t67, t6, t1, 0, 0, 0, 0.2e1 * t57, -0.2e1 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, -t21 * t12 + (-t39 * t7 - t41 * t8) * qJD(5), t20 * t12 + (t39 * t8 - t41 * t7) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, (-t18 * t39 - t19 * t41) * qJD(5) - t21 * qJD(4), (-t18 * t41 + t19 * t39) * qJD(5) + t20 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
