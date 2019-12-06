% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:57
% EndTime: 2019-12-05 17:59:59
% DurationCPUTime: 0.29s
% Computational Cost: add. (419->62), mult. (842->104), div. (0->0), fcn. (714->4), ass. (0->42)
t25 = sin(qJ(4));
t26 = sin(qJ(3));
t27 = cos(qJ(4));
t28 = cos(qJ(3));
t13 = t25 * t28 + t27 * t26;
t42 = t27 * t28;
t14 = -t25 * t26 + t42;
t44 = t27 * pkin(3);
t22 = pkin(4) + t44;
t51 = qJD(3) + qJD(4);
t8 = t51 * t13;
t40 = t26 * qJD(3);
t41 = qJD(4) * t25;
t9 = -t25 * t40 - t26 * t41 + t51 * t42;
t52 = -t22 * t8 + (t25 * t9 + (t13 * t27 - t14 * t25) * qJD(4)) * pkin(3);
t29 = -pkin(1) - pkin(6);
t43 = pkin(7) - t29;
t11 = t43 * t40;
t16 = t43 * t28;
t12 = qJD(3) * t16;
t15 = t43 * t26;
t31 = -t25 * t15 + t27 * t16;
t3 = t31 * qJD(4) - t25 * t11 + t27 * t12;
t49 = 2 * qJD(2);
t48 = pkin(4) * t8;
t47 = t14 * t8;
t20 = t26 * pkin(3) + qJ(2);
t39 = t28 * qJD(3);
t17 = pkin(3) * t39 + qJD(2);
t38 = qJ(2) * qJD(3);
t37 = pkin(3) * t41;
t36 = qJD(4) * t44;
t35 = t13 * t9 - t47;
t32 = t27 * t15 + t25 * t16;
t1 = -t9 * qJ(5) - t13 * qJD(5) - t3;
t4 = t32 * qJD(4) + t27 * t11 + t25 * t12;
t2 = t8 * qJ(5) - t14 * qJD(5) + t4;
t5 = -t14 * qJ(5) - t31;
t6 = -t13 * qJ(5) - t32;
t30 = t1 * t13 + t2 * t14 - t5 * t8 + t6 * t9;
t7 = t9 * pkin(4) + t17;
t10 = [0, 0, 0, 0, t49, qJ(2) * t49, -0.2e1 * t26 * t39, 0.2e1 * (t26 ^ 2 - t28 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t26 + 0.2e1 * t28 * t38, 0.2e1 * qJD(2) * t28 - 0.2e1 * t26 * t38, -0.2e1 * t47, 0.2e1 * t8 * t13 - 0.2e1 * t14 * t9, 0, 0, 0, 0.2e1 * t17 * t13 + 0.2e1 * t20 * t9, 0.2e1 * t17 * t14 - 0.2e1 * t20 * t8, -0.2e1 * t30, 0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * (t13 * pkin(4) + t20) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t35, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t35; 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t39, 0, -t29 * t40, -t29 * t39, 0, 0, -t8, -t9, 0, t4, t3, -t52, t2 * t22 + (t1 * t25 + (-t25 * t5 + t27 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t39, 0, 0, 0, 0, 0, -t8, -t9, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t37, -0.2e1 * t36, 0, 0.2e1 * (-t22 + t44) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, 0, t4, t3, t48, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t36, 0, -pkin(4) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
