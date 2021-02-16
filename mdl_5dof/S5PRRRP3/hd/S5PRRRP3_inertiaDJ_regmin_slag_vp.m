% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:27
% EndTime: 2021-01-15 16:23:30
% DurationCPUTime: 0.31s
% Computational Cost: add. (437->67), mult. (1071->119), div. (0->0), fcn. (912->4), ass. (0->41)
t51 = qJD(3) + qJD(4);
t50 = pkin(6) + pkin(7);
t30 = sin(qJ(4));
t31 = sin(qJ(3));
t32 = cos(qJ(3));
t48 = cos(qJ(4));
t20 = t30 * t32 + t48 * t31;
t12 = t51 * t20;
t49 = t12 * pkin(4);
t37 = t48 * qJD(4);
t38 = t48 * qJD(3);
t46 = t30 * t31;
t11 = t51 * t46 + (-t38 - t37) * t32;
t47 = t20 * t11;
t45 = qJD(4) * t30;
t44 = t31 * qJD(3);
t43 = t32 * qJD(3);
t42 = -0.2e1 * pkin(2) * qJD(3);
t41 = pkin(3) * t44;
t40 = pkin(3) * t45;
t39 = t48 * pkin(3);
t29 = -t32 * pkin(3) - pkin(2);
t36 = qJD(3) * t30 * t50;
t35 = pkin(3) * t37;
t34 = t50 * t38;
t21 = t50 * t31;
t22 = t50 * t32;
t33 = t30 * t21 - t48 * t22;
t3 = t21 * t37 + t22 * t45 + t31 * t34 + t32 * t36;
t4 = t33 * qJD(4) + t31 * t36 - t32 * t34;
t28 = t39 + pkin(4);
t26 = -0.2e1 * t35;
t25 = -0.2e1 * t40;
t19 = -t48 * t32 + t46;
t13 = t19 * pkin(4) + t29;
t9 = t41 + t49;
t8 = -t19 * qJ(5) - t33;
t7 = -t20 * qJ(5) - t48 * t21 - t30 * t22;
t2 = t11 * qJ(5) - t20 * qJD(5) + t4;
t1 = t12 * qJ(5) + t19 * qJD(5) + t3;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t19 * t12 - 0.2e1 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20 * t1 - t11 * t8 - t12 * t7 - t19 * t2; 0, 0, 0, 0, 0.2e1 * t31 * t43, 0.2e1 * (-t31 ^ 2 + t32 ^ 2) * qJD(3), 0, 0, 0, t31 * t42, t32 * t42, -0.2e1 * t47, 0.2e1 * t19 * t11 - 0.2e1 * t20 * t12, 0, 0, 0, 0.2e1 * t29 * t12 + 0.2e1 * t19 * t41, -0.2e1 * t29 * t11 + 0.2e1 * t20 * t41, 0.2e1 * t13 * t12 + 0.2e1 * t9 * t19, -0.2e1 * t13 * t11 + 0.2e1 * t9 * t20, 0.2e1 * t1 * t19 + 0.2e1 * t7 * t11 - 0.2e1 * t8 * t12 - 0.2e1 * t2 * t20, -0.2e1 * t8 * t1 + 0.2e1 * t13 * t9 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t43, 0, 0, 0, 0, 0, -t12, t11, -t12, t11, 0, -t12 * t28 + (-t11 * t30 + (t19 * t30 + t48 * t20) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, t43, -t44, 0, -pkin(6) * t43, pkin(6) * t44, 0, 0, -t11, -t12, 0, t4, t3, t2, t1, t28 * t11 + (-t12 * t30 + (-t48 * t19 + t20 * t30) * qJD(4)) * pkin(3), t2 * t28 + (-t1 * t30 + (-t30 * t7 + t48 * t8) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26, t25, t26, 0, 0.2e1 * (t39 - t28) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, -t12, t11, 0, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, t4, t3, t2, t1, t11 * pkin(4), t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t35, -t40, -t35, 0, -pkin(4) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
