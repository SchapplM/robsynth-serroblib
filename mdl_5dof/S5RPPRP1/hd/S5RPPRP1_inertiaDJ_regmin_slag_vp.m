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
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:36:25
% EndTime: 2019-12-05 17:36:26
% DurationCPUTime: 0.23s
% Computational Cost: add. (183->42), mult. (403->85), div. (0->0), fcn. (313->6), ass. (0->33)
t16 = sin(qJ(4));
t14 = sin(pkin(8));
t33 = qJ(5) * t14;
t15 = cos(pkin(8));
t8 = -cos(pkin(7)) * pkin(1) - pkin(2) - t14 * pkin(6) - t15 * pkin(3);
t22 = -t8 + t33;
t11 = sin(pkin(7)) * pkin(1) + qJ(3);
t17 = cos(qJ(4));
t27 = t17 * t15 * t11;
t38 = t22 * t16 - t27;
t28 = qJD(5) * t14;
t34 = t11 * t16;
t29 = qJD(4) * t17;
t31 = qJD(3) * t17;
t35 = t15 * t31 + t8 * t29;
t1 = -t16 * t28 + (-t15 * t34 - t17 * t33) * qJD(4) + t35;
t32 = qJD(3) * t16;
t23 = t15 * t32;
t2 = t38 * qJD(4) - t17 * t28 - t23;
t3 = -t22 * t17 + (-pkin(4) - t34) * t15;
t37 = -t1 * t16 - t2 * t17 + (t16 * t3 + t17 * t38) * qJD(4);
t36 = 0.2e1 * qJD(4);
t30 = qJD(4) * t16;
t26 = t14 * t30;
t25 = t15 * t30;
t24 = t14 * t29;
t21 = t14 * t15 * t36;
t12 = t14 ^ 2;
t20 = 0.2e1 * (t15 ^ 2 + t12) * qJD(3);
t9 = (pkin(4) * t29 + qJD(3)) * t14;
t5 = -t23 + (-t16 * t8 - t27) * qJD(4);
t4 = t11 * t25 - t35;
t6 = [0, 0, 0, 0, 0, 0, t20, t11 * t20, -0.2e1 * t12 * t16 * t29, (t16 ^ 2 - t17 ^ 2) * t12 * t36, t16 * t21, t17 * t21, 0, -0.2e1 * t5 * t15 + 0.2e1 * (t11 * t29 + t32) * t12, -0.2e1 * t4 * t15 + 0.2e1 * (-t11 * t30 + t31) * t12, 0.2e1 * t37 * t14, -0.2e1 * t38 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * (pkin(4) * t16 + t11) * t9 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t15 + (t1 * t17 - t16 * t2 + (t16 * t38 - t17 * t3) * qJD(4)) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t15 * t29, 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t24, 0, t5, t4, pkin(4) * t26, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t26, 0, -pkin(4) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t29, 0, -pkin(4) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
