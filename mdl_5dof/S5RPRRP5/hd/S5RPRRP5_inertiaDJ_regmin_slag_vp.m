% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP5
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
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:57
% EndTime: 2019-12-31 18:40:58
% DurationCPUTime: 0.25s
% Computational Cost: add. (223->46), mult. (531->77), div. (0->0), fcn. (347->6), ass. (0->38)
t28 = sin(qJ(3));
t33 = cos(pkin(8)) * pkin(1) + pkin(2);
t40 = cos(qJ(3));
t42 = pkin(1) * sin(pkin(8));
t45 = t28 * t42 - t40 * t33;
t44 = 2 * qJD(5);
t30 = t28 * t33 + t40 * t42;
t11 = t30 * qJD(3);
t12 = -pkin(3) + t45;
t29 = cos(qJ(4));
t23 = t29 * qJD(4);
t27 = sin(qJ(4));
t43 = t11 * t27 + t12 * t23;
t39 = t27 * qJD(4);
t15 = -pkin(4) * t39 + qJ(5) * t23 + t27 * qJD(5);
t4 = t11 - t15;
t41 = t15 - t4;
t38 = pkin(3) * t39;
t37 = pkin(3) * t23;
t36 = pkin(7) * t39;
t35 = pkin(7) * t23;
t34 = -t11 * t29 + t12 * t39;
t10 = t45 * qJD(3);
t24 = t27 ^ 2;
t25 = t29 ^ 2;
t1 = (-t24 - t25) * t10;
t32 = -t29 * pkin(4) - t27 * qJ(5);
t14 = t32 * qJD(4) + t29 * qJD(5);
t20 = 0.2e1 * t27 * t23;
t19 = -pkin(3) + t32;
t17 = 0.2e1 * (-t24 + t25) * qJD(4);
t16 = t19 * t39;
t13 = pkin(7) + t30;
t6 = t12 + t32;
t5 = t6 * t39;
t3 = -t27 * t10 + t13 * t23;
t2 = t29 * t10 + t13 * t39;
t7 = [0, 0, 0, 0, 0, -0.2e1 * t11, 0.2e1 * t10, t20, t17, 0, 0, 0, 0.2e1 * t34, 0.2e1 * t43, -0.2e1 * t4 * t29 + 0.2e1 * t5, 0.2e1 * t1, -0.2e1 * t23 * t6 - 0.2e1 * t4 * t27, 0.2e1 * t1 * t13 + 0.2e1 * t6 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t11, t10, t20, t17, 0, 0, 0, t34 - t38, -t37 + t43, t41 * t29 + t16 + t5, t1, t41 * t27 + (-t19 - t6) * t23, pkin(7) * t1 - t6 * t15 + t4 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t20, t17, 0, 0, 0, -0.2e1 * t38, -0.2e1 * t37, 0.2e1 * t15 * t29 + 0.2e1 * t16, 0, 0.2e1 * t15 * t27 - 0.2e1 * t19 * t23, -0.2e1 * t19 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t39, 0, -t3, t2, -t3, t14, -t2, (pkin(4) * t27 - qJ(5) * t29) * t10 + t14 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t23, -t39, 0, t23, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t39, 0, -t35, t36, -t35, t14, -t36, t14 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, qJ(5) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
