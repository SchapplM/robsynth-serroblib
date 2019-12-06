% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP2
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
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:44
% EndTime: 2019-12-05 18:01:45
% DurationCPUTime: 0.21s
% Computational Cost: add. (264->62), mult. (580->100), div. (0->0), fcn. (402->6), ass. (0->44)
t24 = cos(pkin(8)) * pkin(1) + pkin(2);
t32 = sin(qJ(3));
t34 = cos(qJ(3));
t48 = pkin(1) * sin(pkin(8));
t50 = -t34 * t24 + t32 * t48;
t35 = t32 * t24 + t34 * t48;
t13 = t35 * qJD(3);
t14 = -pkin(3) + t50;
t33 = cos(qJ(4));
t28 = qJD(4) * t33;
t31 = sin(qJ(4));
t49 = t13 * t31 + t14 * t28;
t47 = t33 * pkin(4);
t12 = t50 * qJD(3);
t46 = t33 * t12;
t44 = -qJ(5) - pkin(7);
t15 = pkin(7) + t35;
t43 = -qJ(5) - t15;
t42 = qJD(4) * t31;
t41 = pkin(3) * t42;
t40 = pkin(3) * t28;
t26 = pkin(4) * t42;
t39 = pkin(4) * t28;
t38 = -t13 * t33 + t14 * t42;
t37 = qJD(4) * t44;
t36 = qJD(4) * t43;
t29 = t33 * qJ(5);
t27 = t33 * qJD(5);
t25 = -pkin(3) - t47;
t22 = 0.2e1 * t31 * t28;
t21 = t33 * pkin(7) + t29;
t20 = t44 * t31;
t18 = 0.2e1 * (-t31 ^ 2 + t33 ^ 2) * qJD(4);
t17 = -t31 * qJD(5) + t33 * t37;
t16 = t31 * t37 + t27;
t11 = t16 * t33;
t10 = t14 - t47;
t6 = t26 + t13;
t5 = t33 * t15 + t29;
t4 = t43 * t31;
t3 = (-qJD(5) + t12) * t31 + t33 * t36;
t2 = t31 * t36 + t27 - t46;
t1 = t2 * t33;
t7 = [0, 0, 0, 0, 0, -0.2e1 * t13, 0.2e1 * t12, t22, t18, 0, 0, 0, 0.2e1 * t38, 0.2e1 * t49, -0.2e1 * t3 * t31 + 0.2e1 * t1 + 0.2e1 * (-t31 * t5 - t33 * t4) * qJD(4), 0.2e1 * t10 * t6 + 0.2e1 * t5 * t2 + 0.2e1 * t4 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * t31 + t3 * t33 + (-t31 * t4 + t33 * t5) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t13, t12, t22, t18, 0, 0, 0, t38 - t41, -t40 + t49, t1 + t11 + (-t17 - t3) * t31 + ((-t20 - t4) * t33 + (-t21 - t5) * t31) * qJD(4), t10 * t26 + t5 * t16 + t4 * t17 + t2 * t21 + t3 * t20 + t6 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * t16 + t33 * t17 + (-t20 * t31 + t21 * t33) * qJD(4); 0, 0, 0, 0, 0, 0, 0, t22, t18, 0, 0, 0, -0.2e1 * t41, -0.2e1 * t40, -0.2e1 * t17 * t31 + 0.2e1 * t11 + 0.2e1 * (-t20 * t33 - t21 * t31) * qJD(4), 0.2e1 * t21 * t16 + 0.2e1 * t20 * t17 + 0.2e1 * t25 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t42, 0, t31 * t12 - t15 * t28, t15 * t42 + t46, -t39, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t28, 0, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t42, 0, -pkin(7) * t28, pkin(7) * t42, -t39, t17 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
