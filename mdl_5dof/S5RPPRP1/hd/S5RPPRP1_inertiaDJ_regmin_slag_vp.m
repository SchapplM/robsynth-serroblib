% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:45
% EndTime: 2022-01-23 09:12:46
% DurationCPUTime: 0.25s
% Computational Cost: add. (229->55), mult. (515->100), div. (0->0), fcn. (403->6), ass. (0->38)
t17 = sin(pkin(8));
t43 = (-qJ(5) - pkin(6)) * t17;
t18 = cos(pkin(8));
t19 = sin(qJ(4));
t20 = cos(qJ(4));
t33 = t17 * qJD(5);
t14 = sin(pkin(7)) * pkin(1) + qJ(3);
t38 = t19 * t14;
t24 = -cos(pkin(7)) * pkin(1) - pkin(2);
t40 = t18 * pkin(3);
t23 = t24 - t40;
t21 = -t17 * pkin(6) + t23;
t34 = qJD(4) * t20;
t36 = qJD(3) * t20;
t39 = -t18 * t36 - t21 * t34;
t1 = t19 * t33 + (t20 * t17 * qJ(5) + t18 * t38) * qJD(4) + t39;
t32 = t19 * qJD(3);
t29 = t18 * t32;
t31 = t20 * t18 * t14;
t2 = -t29 - t20 * t33 + (-t31 + (-t23 - t43) * t19) * qJD(4);
t22 = t24 + t43;
t3 = t22 * t20 + (-t20 * pkin(3) - pkin(4) - t38) * t18;
t6 = t31 + (t22 - t40) * t19;
t42 = t1 * t19 - t2 * t20 + (t19 * t3 - t20 * t6) * qJD(4);
t41 = 0.2e1 * qJD(4);
t35 = qJD(4) * t19;
t11 = t17 * t35;
t12 = t18 * t35;
t30 = t17 * t34;
t28 = t17 * t18 * t41;
t15 = t17 ^ 2;
t27 = 0.2e1 * (t18 ^ 2 + t15) * qJD(3);
t13 = t18 * t34;
t9 = (pkin(4) * t34 + qJD(3)) * t17;
t8 = (pkin(4) * t19 + t14) * t17;
t5 = -t29 + (-t19 * t21 - t31) * qJD(4);
t4 = t14 * t12 + t39;
t7 = [0, 0, 0, 0, 0, t27, t14 * t27, -0.2e1 * t15 * t19 * t34, (t19 ^ 2 - t20 ^ 2) * t15 * t41, t19 * t28, t20 * t28, 0, -0.2e1 * t5 * t18 + 0.2e1 * (t14 * t34 + t32) * t15, -0.2e1 * t4 * t18 + 0.2e1 * (-t14 * t35 + t36) * t15, -0.2e1 * t2 * t18 + 0.2e1 * (t19 * t9 + t8 * t34) * t17, -0.2e1 * t1 * t18 + 0.2e1 * (t20 * t9 - t8 * t35) * t17, 0.2e1 * t42 * t17, -0.2e1 * t6 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t18 + (-t1 * t20 - t19 * t2 + (-t19 * t6 - t20 * t3) * qJD(4)) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13, t12, t13, 0, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t30, 0, t5, t4, t2, t1, pkin(4) * t11, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t11, -t30, t11, 0, -pkin(4) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t34, -t35, -t34, 0, -pkin(4) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t11, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
