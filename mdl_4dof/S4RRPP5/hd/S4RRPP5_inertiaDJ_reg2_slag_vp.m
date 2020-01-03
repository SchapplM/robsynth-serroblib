% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_inertiaDJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:37
% EndTime: 2019-12-31 17:00:38
% DurationCPUTime: 0.23s
% Computational Cost: add. (81->38), mult. (237->77), div. (0->0), fcn. (127->2), ass. (0->31)
t20 = cos(qJ(2));
t16 = t20 * qJD(3);
t19 = sin(qJ(2));
t29 = t19 * qJ(3);
t23 = -t20 * pkin(2) - t29;
t33 = t23 * qJD(2) + t16;
t18 = -pkin(2) - qJ(4);
t32 = t18 * t20 - t29;
t31 = pkin(3) + pkin(5);
t15 = t19 * qJD(2);
t17 = t20 * qJD(2);
t28 = qJ(3) * qJD(3);
t27 = -0.2e1 * pkin(1) * qJD(2);
t26 = pkin(5) * t15;
t25 = t19 * t17;
t24 = pkin(2) * t15 - t19 * qJD(3);
t21 = 0.2e1 * qJD(3);
t13 = pkin(5) * t17;
t12 = -0.2e1 * t25;
t11 = 0.2e1 * t25;
t10 = t31 * t20;
t9 = t31 * t19;
t8 = (-t19 ^ 2 + t20 ^ 2) * qJD(2);
t7 = -pkin(1) + t23;
t6 = pkin(3) * t17 + t13;
t5 = t31 * t15;
t4 = 0.2e1 * t8;
t3 = -pkin(1) + t32;
t2 = -qJ(3) * t17 + t24;
t1 = -t20 * qJD(4) + (-qJ(3) * t20 + qJ(4) * t19) * qJD(2) + t24;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t4, 0, t12, 0, 0, t19 * t27, t20 * t27, 0, 0, 0, 0, 0, t11, t4, t12, 0, -0.2e1 * t7 * t15 + 0.2e1 * t2 * t20, -0.2e1 * t7 * t17 - 0.2e1 * t2 * t19, 0.2e1 * t7 * t2, 0, 0, 0, t12, -0.2e1 * t8, t11, 0.2e1 * t6 * t19 - 0.2e1 * t5 * t20 + 0.2e1 * (-t10 * t19 + t20 * t9) * qJD(2), -0.2e1 * t1 * t19 - 0.2e1 * t3 * t17, -0.2e1 * t1 * t20 + 0.2e1 * t3 * t15, 0.2e1 * t3 * t1 - 0.2e1 * t10 * t5 + 0.2e1 * t9 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t15, 0, -t13, t26, 0, 0, 0, -t17, t15, 0, 0, 0, t33, t13, -t26, t33 * pkin(5), 0, t15, t17, 0, 0, 0, t32 * qJD(2) - qJD(4) * t19 + t16, -t5, -t6, -t5 * qJ(3) + t10 * qJD(3) - t9 * qJD(4) + t6 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0.2e1 * t28, 0, 0, 0, 0, 0, 0, 0, t21, 0.2e1 * qJD(4), -0.2e1 * t18 * qJD(4) + 0.2e1 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, t13, 0, 0, 0, 0, 0, 0, t17, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
