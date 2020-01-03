% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:06
% EndTime: 2019-12-31 16:35:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (106->33), mult. (321->74), div. (0->0), fcn. (276->6), ass. (0->33)
t32 = qJD(3) + qJD(4);
t14 = sin(qJ(4));
t15 = sin(qJ(3));
t17 = cos(qJ(4));
t18 = cos(qJ(3));
t7 = t14 * t15 - t17 * t18;
t5 = t32 * t7;
t31 = pkin(5) + pkin(6);
t30 = pkin(3) * qJD(4);
t16 = sin(qJ(2));
t29 = qJD(2) * t16;
t19 = cos(qJ(2));
t28 = qJD(2) * t19;
t27 = qJD(3) * t15;
t26 = qJD(3) * t18;
t25 = qJD(3) * t19;
t24 = -0.2e1 * pkin(2) * qJD(3);
t23 = pkin(3) * t27;
t22 = t14 * t30;
t21 = t17 * t30;
t20 = qJD(3) * t31;
t8 = t14 * t18 + t17 * t15;
t6 = t32 * t8;
t13 = -t18 * pkin(3) - pkin(2);
t12 = t31 * t18;
t11 = t31 * t15;
t10 = t18 * t20;
t9 = t15 * t20;
t4 = -t17 * t10 + t14 * t9 + (t11 * t14 - t12 * t17) * qJD(4);
t3 = t14 * t10 + t17 * t9 + (t11 * t17 + t12 * t14) * qJD(4);
t2 = t16 * t5 - t8 * t28;
t1 = t6 * t16 + t7 * t28;
t33 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t29, -t28, 0, 0, 0, 0, 0, -t15 * t25 - t18 * t29, t15 * t29 - t18 * t25, 0, 0, 0, 0, 0, -t19 * t6 + t7 * t29, t19 * t5 + t8 * t29; 0, 0, 0, 0, 0.2e1 * t15 * t26, 0.2e1 * (-t15 ^ 2 + t18 ^ 2) * qJD(3), 0, 0, 0, t15 * t24, t18 * t24, -0.2e1 * t8 * t5, 0.2e1 * t5 * t7 - 0.2e1 * t8 * t6, 0, 0, 0, 0.2e1 * t13 * t6 + 0.2e1 * t7 * t23, -0.2e1 * t13 * t5 + 0.2e1 * t8 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t28 - t16 * t26, t16 * t27 - t18 * t28, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, t26, -t27, 0, -pkin(5) * t26, pkin(5) * t27, 0, 0, -t5, -t6, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t22, -0.2e1 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t33;
