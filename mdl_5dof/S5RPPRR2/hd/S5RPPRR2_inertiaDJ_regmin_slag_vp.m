% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:56
% EndTime: 2019-12-05 17:39:58
% DurationCPUTime: 0.30s
% Computational Cost: add. (316->43), mult. (695->88), div. (0->0), fcn. (674->6), ass. (0->38)
t30 = sin(pkin(8));
t31 = cos(pkin(8));
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t53 = t34 * t30 - t31 * t36;
t23 = (t30 ^ 2 + t31 ^ 2) * qJD(3);
t19 = t30 * t36 + t34 * t31;
t17 = t19 * qJD(4);
t18 = t53 * qJD(4);
t33 = sin(qJ(5));
t35 = cos(qJ(5));
t41 = t19 * t33 + t35 * t53;
t3 = t41 * qJD(5) + t33 * t17 + t35 * t18;
t32 = -pkin(1) - qJ(3);
t50 = -pkin(6) + t32;
t21 = t50 * t30;
t22 = t50 * t31;
t39 = t34 * t21 - t22 * t36;
t52 = t19 * qJD(3) + t39 * qJD(4);
t51 = 2 * qJD(2);
t46 = pkin(4) * qJD(5);
t26 = t30 * pkin(3) + qJ(2);
t45 = qJ(2) * qJD(2);
t44 = t33 * t46;
t43 = t35 * t46;
t11 = t19 * t35 - t33 * t53;
t40 = -t21 * t36 - t34 * t22;
t38 = -t11 * qJD(5) - t35 * t17 + t33 * t18;
t37 = t53 * qJD(3) + t40 * qJD(4);
t14 = -pkin(4) * t18 + qJD(2);
t13 = pkin(4) * t19 + t26;
t10 = -t19 * pkin(7) - t40;
t9 = pkin(7) * t53 - t39;
t8 = t17 * pkin(7) + t37;
t7 = t18 * pkin(7) - t52;
t2 = -t33 * t7 + t35 * t8 + (-t10 * t35 - t33 * t9) * qJD(5);
t1 = -t33 * t8 - t35 * t7 + (t10 * t33 - t35 * t9) * qJD(5);
t4 = [0, 0, 0, 0, t51, 0.2e1 * t45, t30 * t51, t31 * t51, 0.2e1 * t23, -0.2e1 * t32 * t23 + 0.2e1 * t45, 0.2e1 * t53 * t17, 0.2e1 * t17 * t19 - 0.2e1 * t18 * t53, 0, 0, 0, 0.2e1 * qJD(2) * t19 - 0.2e1 * t18 * t26, -0.2e1 * qJD(2) * t53 - 0.2e1 * t17 * t26, -0.2e1 * t41 * t38, -0.2e1 * t11 * t38 - 0.2e1 * t3 * t41, 0, 0, 0, 0.2e1 * t11 * t14 - 0.2e1 * t13 * t3, 0.2e1 * t13 * t38 - 0.2e1 * t14 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, -t18, -t17, 0, 0, 0, 0, 0, -t3, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t18, 0, t37, t52, 0, 0, t38, t3, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t18, 0, 0, 0, 0, 0, t38, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t44, -0.2e1 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t3, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
