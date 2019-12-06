% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP1
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
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:12
% EndTime: 2019-12-05 16:40:14
% DurationCPUTime: 0.21s
% Computational Cost: add. (160->60), mult. (395->97), div. (0->0), fcn. (243->4), ass. (0->41)
t28 = cos(qJ(3));
t43 = t28 * pkin(2);
t42 = -qJ(5) - pkin(7);
t19 = -pkin(3) - t43;
t27 = cos(qJ(4));
t23 = qJD(4) * t27;
t25 = sin(qJ(4));
t26 = sin(qJ(3));
t40 = pkin(2) * qJD(3);
t35 = t26 * t40;
t41 = t19 * t23 + t25 * t35;
t18 = t26 * pkin(2) + pkin(7);
t39 = -qJ(5) - t18;
t38 = qJD(4) * t25;
t37 = pkin(3) * t38;
t36 = pkin(3) * t23;
t34 = t28 * t40;
t21 = pkin(4) * t38;
t33 = pkin(4) * t23;
t20 = -t27 * pkin(4) - pkin(3);
t32 = qJD(4) * t42;
t31 = qJD(4) * t39;
t30 = t27 * t34;
t29 = t19 * t38 - t27 * t35;
t24 = t27 * qJ(5);
t22 = t27 * qJD(5);
t17 = 0.2e1 * t25 * t23;
t15 = t27 * pkin(7) + t24;
t14 = t42 * t25;
t13 = t20 - t43;
t10 = t21 + t35;
t9 = 0.2e1 * (-t25 ^ 2 + t27 ^ 2) * qJD(4);
t8 = t27 * t18 + t24;
t7 = t39 * t25;
t6 = -t25 * qJD(5) + t27 * t32;
t5 = t25 * t32 + t22;
t4 = t5 * t27;
t3 = (-qJD(5) - t34) * t25 + t27 * t31;
t2 = t25 * t31 + t22 + t30;
t1 = t2 * t27;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * t2 + t27 * t3 + (-t25 * t7 + t27 * t8) * qJD(4); 0, 0, 0, 0, 0, -0.2e1 * t35, -0.2e1 * t34, t17, t9, 0, 0, 0, 0.2e1 * t29, 0.2e1 * t41, -0.2e1 * t3 * t25 + 0.2e1 * t1 + 0.2e1 * (-t25 * t8 - t27 * t7) * qJD(4), 0.2e1 * t13 * t10 + 0.2e1 * t8 * t2 + 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * t5 + t27 * t6 + (-t14 * t25 + t15 * t27) * qJD(4); 0, 0, 0, 0, 0, -t35, -t34, t17, t9, 0, 0, 0, t29 - t37, -t36 + t41, t1 + t4 + (-t3 - t6) * t25 + ((-t14 - t7) * t27 + (-t15 - t8) * t25) * qJD(4), t10 * t20 + t13 * t21 + t3 * t14 + t2 * t15 + t8 * t5 + t7 * t6; 0, 0, 0, 0, 0, 0, 0, t17, t9, 0, 0, 0, -0.2e1 * t37, -0.2e1 * t36, -0.2e1 * t6 * t25 + 0.2e1 * t4 + 0.2e1 * (-t14 * t27 - t15 * t25) * qJD(4), 0.2e1 * t14 * t6 + 0.2e1 * t15 * t5 + 0.2e1 * t20 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t23, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t38, 0, -t18 * t23 - t25 * t34, t18 * t38 - t30, -t33, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t38, 0, -pkin(7) * t23, pkin(7) * t38, -t33, t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
