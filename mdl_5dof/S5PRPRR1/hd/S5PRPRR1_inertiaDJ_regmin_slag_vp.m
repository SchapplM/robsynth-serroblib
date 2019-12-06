% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:04
% EndTime: 2019-12-05 15:43:05
% DurationCPUTime: 0.21s
% Computational Cost: add. (270->46), mult. (674->82), div. (0->0), fcn. (662->6), ass. (0->39)
t27 = sin(pkin(9));
t28 = cos(pkin(9));
t30 = sin(qJ(4));
t31 = cos(qJ(4));
t33 = t30 * t27 - t31 * t28;
t45 = pkin(6) + qJ(3);
t19 = t45 * t27;
t20 = t45 * t28;
t35 = -t31 * t19 - t30 * t20;
t48 = t33 * qJD(3) - t35 * qJD(4);
t42 = qJD(4) * t31;
t43 = qJD(4) * t30;
t16 = t27 * t43 - t28 * t42;
t47 = -0.2e1 * t16;
t46 = cos(qJ(5));
t44 = t27 * t42 + t28 * t43;
t29 = sin(qJ(5));
t41 = qJD(5) * t29;
t40 = pkin(4) * t41;
t24 = -t28 * pkin(3) - pkin(2);
t39 = pkin(4) * t44;
t38 = qJD(5) * t46;
t37 = pkin(4) * t38;
t36 = 0.2e1 * (t27 ^ 2 + t28 ^ 2) * qJD(3);
t34 = t30 * t19 - t31 * t20;
t18 = t31 * t27 + t30 * t28;
t10 = t46 * t18 - t29 * t33;
t32 = -t18 * qJD(3) + t34 * qJD(4);
t11 = pkin(4) * t33 + t24;
t9 = t29 * t18 + t33 * t46;
t8 = -pkin(7) * t33 - t34;
t7 = -t18 * pkin(7) + t35;
t6 = t16 * pkin(7) + t32;
t5 = -t44 * pkin(7) - t48;
t4 = t10 * qJD(5) - t29 * t16 + t46 * t44;
t3 = t46 * t16 + t18 * t41 + t29 * t44 + t33 * t38;
t2 = t46 * t6 - t29 * t5 + (-t29 * t7 - t46 * t8) * qJD(5);
t1 = -t46 * t5 - t29 * t6 + (t29 * t8 - t46 * t7) * qJD(5);
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t36, qJ(3) * t36, t18 * t47, 0.2e1 * t16 * t33 - 0.2e1 * t18 * t44, 0, 0, 0, 0.2e1 * t24 * t44, t24 * t47, -0.2e1 * t10 * t3, -0.2e1 * t10 * t4 + 0.2e1 * t3 * t9, 0, 0, 0, 0.2e1 * t11 * t4 + 0.2e1 * t9 * t39, 0.2e1 * t10 * t39 - 0.2e1 * t11 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t16, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t16, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t44, 0, t32, t48, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t40, -0.2e1 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
