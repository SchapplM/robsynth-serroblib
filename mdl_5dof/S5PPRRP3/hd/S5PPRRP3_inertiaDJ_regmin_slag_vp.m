% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:16
% EndTime: 2019-12-05 15:11:17
% DurationCPUTime: 0.20s
% Computational Cost: add. (117->49), mult. (409->93), div. (0->0), fcn. (349->6), ass. (0->41)
t41 = 2 * qJD(5);
t17 = sin(pkin(8));
t22 = cos(qJ(3));
t40 = t17 * t22;
t19 = sin(qJ(4));
t15 = t19 ^ 2;
t21 = cos(qJ(4));
t16 = t21 ^ 2;
t39 = t15 + t16;
t20 = sin(qJ(3));
t38 = qJD(3) * t20;
t37 = qJD(3) * t22;
t36 = qJD(4) * t19;
t14 = qJD(4) * t21;
t35 = qJD(4) * t22;
t34 = -0.2e1 * pkin(3) * qJD(4);
t33 = pkin(6) * t36;
t32 = pkin(6) * t14;
t31 = t17 * t38;
t30 = t21 * t38;
t29 = t20 * t37;
t28 = t21 * t35;
t27 = t39 * t22;
t26 = -t21 * pkin(4) - t19 * qJ(5);
t25 = pkin(4) * t19 - qJ(5) * t21;
t18 = cos(pkin(8));
t6 = t18 * t21 + t19 * t40;
t8 = t20 * t36 - t21 * t37;
t9 = t20 * t14 + t19 * t37;
t24 = t26 * qJD(4) + t21 * qJD(5);
t1 = -t17 * t28 + (qJD(4) * t18 + t31) * t19;
t2 = qJD(4) * t6 + t17 * t30;
t7 = -t18 * t19 + t21 * t40;
t23 = -t1 * t19 - t2 * t21 + (-t19 * t7 + t21 * t6) * qJD(4);
t12 = -pkin(3) + t26;
t11 = -t19 * t35 - t30;
t10 = t19 * t38 - t28;
t5 = t25 * qJD(4) - t19 * qJD(5);
t4 = t8 * t17;
t3 = t9 * t17;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t17 ^ 2 * t29 - 0.2e1 * t6 * t1 - 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t19 * t6 + t21 * t7 - t40) * t37 + (t23 + t31) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t39) * t29; 0, 0, 0, -t17 * t37, t31, 0, 0, 0, 0, 0, t4, t3, t4, t23, -t3, (t12 * t37 + t20 * t5) * t17 + t23 * pkin(6); 0, 0, 0, -t38, -t37, 0, 0, 0, 0, 0, t11, t10, t11, qJD(3) * t27, -t10, -t22 * t5 + (pkin(6) * t27 + t12 * t20) * qJD(3); 0, 0, 0, 0, 0, 0.2e1 * t19 * t14, 0.2e1 * (-t15 + t16) * qJD(4), 0, 0, 0, t19 * t34, t21 * t34, 0.2e1 * t12 * t36 - 0.2e1 * t5 * t21, 0, -0.2e1 * t12 * t14 - 0.2e1 * t5 * t19, 0.2e1 * t12 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, t1 * pkin(4) - t2 * qJ(5) + t7 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8, -t9, 0, -t8, t24 * t20 - t25 * t37; 0, 0, 0, 0, 0, 0, 0, t14, -t36, 0, -t32, t33, -t32, t24, -t33, t24 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, qJ(5) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
