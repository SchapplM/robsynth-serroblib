% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:59
% EndTime: 2019-12-05 15:17:01
% DurationCPUTime: 0.31s
% Computational Cost: add. (159->47), mult. (495->98), div. (0->0), fcn. (460->8), ass. (0->44)
t22 = sin(qJ(5));
t23 = sin(qJ(4));
t25 = cos(qJ(5));
t26 = cos(qJ(4));
t14 = t22 * t26 + t25 * t23;
t24 = sin(qJ(3));
t27 = cos(qJ(3));
t39 = t27 * qJD(3);
t13 = t22 * t23 - t25 * t26;
t47 = qJD(4) + qJD(5);
t7 = t47 * t13;
t49 = t14 * t39 - t24 * t7;
t8 = t47 * t14;
t3 = t13 * t39 + t24 * t8;
t46 = pkin(6) + pkin(7);
t20 = sin(pkin(9));
t44 = t20 * t27;
t43 = pkin(4) * qJD(5);
t42 = t23 * qJD(4);
t41 = t24 * qJD(3);
t40 = t26 * qJD(4);
t38 = -0.2e1 * pkin(3) * qJD(4);
t37 = pkin(4) * t42;
t36 = t22 * t43;
t35 = t25 * t43;
t32 = t20 * t41;
t31 = qJD(4) * t46;
t21 = cos(pkin(9));
t12 = -t21 * t23 + t26 * t44;
t30 = t21 * t26 + t23 * t44;
t29 = t24 * t42 - t26 * t39;
t28 = t23 * t39 + t24 * t40;
t19 = -t26 * pkin(4) - pkin(3);
t18 = t46 * t26;
t17 = t46 * t23;
t16 = t26 * t31;
t15 = t23 * t31;
t10 = t30 * qJD(4) + t26 * t32;
t9 = -t12 * qJD(4) + t23 * t32;
t6 = t22 * t15 - t25 * t16 + (t17 * t22 - t18 * t25) * qJD(5);
t5 = t25 * t15 + t22 * t16 + (t17 * t25 + t18 * t22) * qJD(5);
t2 = t22 * t10 + t25 * t9 + (-t12 * t25 + t22 * t30) * qJD(5);
t1 = t25 * t10 - t22 * t9 + (t12 * t22 + t25 * t30) * qJD(5);
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t20 * t39, t32, 0, 0, 0, 0, 0, t29 * t20, t28 * t20, 0, 0, 0, 0, 0, t3 * t20, t49 * t20; 0, 0, 0, -t41, -t39, 0, 0, 0, 0, 0, -t26 * t41 - t27 * t42, t23 * t41 - t27 * t40, 0, 0, 0, 0, 0, t13 * t41 - t27 * t8, t14 * t41 + t27 * t7; 0, 0, 0, 0, 0, 0.2e1 * t23 * t40, 0.2e1 * (-t23 ^ 2 + t26 ^ 2) * qJD(4), 0, 0, 0, t23 * t38, t26 * t38, -0.2e1 * t14 * t7, 0.2e1 * t7 * t13 - 0.2e1 * t14 * t8, 0, 0, 0, 0.2e1 * t13 * t37 + 0.2e1 * t19 * t8, 0.2e1 * t14 * t37 - 0.2e1 * t19 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t29, 0, 0, 0, 0, 0, -t49, t3; 0, 0, 0, 0, 0, 0, 0, t40, -t42, 0, -pkin(6) * t40, pkin(6) * t42, 0, 0, -t7, -t8, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t36, -0.2e1 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
