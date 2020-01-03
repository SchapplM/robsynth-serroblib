% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:57
% EndTime: 2019-12-31 19:29:59
% DurationCPUTime: 0.31s
% Computational Cost: add. (443->73), mult. (1037->135), div. (0->0), fcn. (901->6), ass. (0->53)
t67 = 2 * qJD(4);
t66 = -pkin(3) - pkin(4);
t65 = -qJ(3) - pkin(6);
t50 = sin(qJ(2));
t52 = cos(qJ(2));
t57 = qJD(2) * t65;
t31 = t52 * qJD(3) + t50 * t57;
t47 = sin(pkin(8));
t48 = cos(pkin(8));
t54 = -t50 * qJD(3) + t52 * t57;
t14 = t48 * t31 + t47 * t54;
t39 = t65 * t50;
t40 = t65 * t52;
t23 = t47 * t39 - t48 * t40;
t49 = sin(qJ(5));
t64 = qJD(5) * t49;
t51 = cos(qJ(5));
t63 = qJD(5) * t51;
t62 = t50 * qJD(2);
t61 = t52 * qJD(2);
t60 = -0.2e1 * pkin(1) * qJD(2);
t46 = pkin(2) * t62;
t59 = -t52 * pkin(2) - pkin(1);
t45 = -t48 * pkin(2) - pkin(3);
t13 = t47 * t31 - t48 * t54;
t22 = -t48 * t39 - t47 * t40;
t58 = t22 * t13 + t23 * t14;
t34 = t47 * t50 - t48 * t52;
t35 = t47 * t52 + t48 * t50;
t19 = t49 * t34 + t51 * t35;
t56 = t35 * qJ(4) - t59;
t33 = -t47 * t62 + t48 * t61;
t55 = t33 * qJ(4) + t35 * qJD(4) - t46;
t32 = t35 * qJD(2);
t53 = 0.2e1 * t13 * t35 - 0.2e1 * t14 * t34 + 0.2e1 * t22 * t33 - 0.2e1 * t23 * t32;
t43 = t47 * pkin(2) + qJ(4);
t42 = -pkin(4) + t45;
t21 = t49 * qJD(4) + (t42 * t49 + t43 * t51) * qJD(5);
t20 = -t51 * qJD(4) + (-t42 * t51 + t43 * t49) * qJD(5);
t18 = -t51 * t34 + t49 * t35;
t17 = t34 * pkin(3) - t56;
t16 = t34 * pkin(7) + t23;
t15 = -t35 * pkin(7) + t22;
t11 = t66 * t34 + t56;
t9 = t32 * pkin(3) - t55;
t8 = t32 * pkin(7) + t14;
t7 = -t33 * pkin(7) + t13;
t5 = t66 * t32 + t55;
t4 = t19 * qJD(5) - t51 * t32 + t49 * t33;
t3 = -t49 * t32 - t51 * t33 - t34 * t63 + t35 * t64;
t2 = t49 * t8 - t51 * t7 + (t15 * t49 + t16 * t51) * qJD(5);
t1 = -t49 * t7 - t51 * t8 + (-t15 * t51 + t16 * t49) * qJD(5);
t6 = [0, 0, 0, 0.2e1 * t50 * t61, 0.2e1 * (-t50 ^ 2 + t52 ^ 2) * qJD(2), 0, 0, 0, t50 * t60, t52 * t60, t53, 0.2e1 * t59 * t46 + 0.2e1 * t58, 0.2e1 * t17 * t32 + 0.2e1 * t9 * t34, t53, -0.2e1 * t17 * t33 - 0.2e1 * t9 * t35, 0.2e1 * t17 * t9 + 0.2e1 * t58, -0.2e1 * t19 * t3, 0.2e1 * t3 * t18 - 0.2e1 * t19 * t4, 0, 0, 0, 0.2e1 * t11 * t4 + 0.2e1 * t5 * t18, -0.2e1 * t11 * t3 + 0.2e1 * t5 * t19; 0, 0, 0, 0, 0, t61, -t62, 0, -pkin(6) * t61, pkin(6) * t62, (-t32 * t47 - t33 * t48) * pkin(2), (-t13 * t48 + t14 * t47) * pkin(2), -t13, -qJD(4) * t34 - t43 * t32 + t45 * t33, t14, t23 * qJD(4) + t13 * t45 + t14 * t43, 0, 0, t3, t4, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t43 * t67, 0, 0, 0, 0, 0, 0.2e1 * t21, -0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t32, 0, -t33, t9, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
