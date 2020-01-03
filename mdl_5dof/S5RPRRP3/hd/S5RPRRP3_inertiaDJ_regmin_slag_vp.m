% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:47:43
% EndTime: 2020-01-03 11:47:45
% DurationCPUTime: 0.41s
% Computational Cost: add. (416->61), mult. (910->112), div. (0->0), fcn. (787->6), ass. (0->41)
t48 = qJD(3) + qJD(4);
t28 = sin(qJ(3));
t22 = sin(pkin(8)) * pkin(1) + pkin(6);
t46 = pkin(7) + t22;
t35 = qJD(3) * t46;
t13 = t28 * t35;
t29 = cos(qJ(3));
t14 = t29 * t35;
t27 = sin(qJ(4));
t15 = t46 * t28;
t16 = t46 * t29;
t45 = cos(qJ(4));
t31 = -t45 * t15 - t27 * t16;
t3 = -t31 * qJD(4) + t45 * t13 + t27 * t14;
t18 = t27 * t29 + t45 * t28;
t11 = t48 * t18;
t47 = t11 * pkin(4);
t34 = t45 * qJD(4);
t36 = t45 * t29;
t43 = t27 * t28;
t10 = -qJD(3) * t36 - t29 * t34 + t48 * t43;
t44 = t18 * t10;
t42 = t28 * qJD(3);
t41 = t29 * qJD(3);
t40 = 0.2e1 * t41;
t39 = pkin(3) * t42;
t38 = qJD(4) * t27 * pkin(3);
t37 = t45 * pkin(3);
t23 = -cos(pkin(8)) * pkin(1) - pkin(2);
t33 = pkin(3) * t34;
t19 = -t29 * pkin(3) + t23;
t30 = t27 * t15 - t45 * t16;
t4 = t30 * qJD(4) + t27 * t13 - t45 * t14;
t25 = t37 + pkin(4);
t17 = -t36 + t43;
t9 = t39 + t47;
t6 = -t17 * qJ(5) - t30;
t5 = -t18 * qJ(5) + t31;
t2 = t10 * qJ(5) - t18 * qJD(5) + t4;
t1 = -t11 * qJ(5) - t17 * qJD(5) - t3;
t7 = [0, 0, 0, 0, t28 * t40, 0.2e1 * (-t28 ^ 2 + t29 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t23 * t42, t23 * t40, -0.2e1 * t44, 0.2e1 * t17 * t10 - 0.2e1 * t18 * t11, 0, 0, 0, 0.2e1 * t19 * t11 + 0.2e1 * t17 * t39, -0.2e1 * t19 * t10 + 0.2e1 * t18 * t39, -0.2e1 * t1 * t17 + 0.2e1 * t5 * t10 - 0.2e1 * t6 * t11 - 0.2e1 * t2 * t18, 0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * (t17 * pkin(4) + t19) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t18 - t6 * t10 - t5 * t11 - t2 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t17 * t11 - 0.2e1 * t44; 0, 0, 0, 0, 0, 0, t41, -t42, 0, -t22 * t41, t22 * t42, 0, 0, -t10, -t11, 0, t4, t3, t25 * t10 + (-t11 * t27 + (-t45 * t17 + t18 * t27) * qJD(4)) * pkin(3), t2 * t25 + (t1 * t27 + (-t27 * t5 + t45 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t41, 0, 0, 0, 0, 0, -t11, t10, 0, -t11 * t25 + (-t10 * t27 + (t17 * t27 + t45 * t18) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t38, -0.2e1 * t33, 0, 0.2e1 * (t37 - t25) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, t4, t3, pkin(4) * t10, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t33, 0, -pkin(4) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
