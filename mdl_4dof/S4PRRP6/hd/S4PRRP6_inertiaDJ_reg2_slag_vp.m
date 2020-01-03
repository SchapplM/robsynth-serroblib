% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRP6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:43
% EndTime: 2019-12-31 16:30:43
% DurationCPUTime: 0.21s
% Computational Cost: add. (54->28), mult. (212->56), div. (0->0), fcn. (137->4), ass. (0->32)
t20 = cos(qJ(2));
t29 = t20 * qJD(2);
t17 = sin(qJ(3));
t15 = t17 ^ 2;
t19 = cos(qJ(3));
t16 = t19 ^ 2;
t32 = t15 + t16;
t4 = t32 * t29;
t34 = 2 * qJD(4);
t33 = pkin(5) * t4;
t31 = t17 * qJD(3);
t18 = sin(qJ(2));
t30 = t18 * qJD(2);
t14 = t19 * qJD(3);
t28 = -0.2e1 * pkin(2) * qJD(3);
t26 = pkin(5) * t31;
t25 = pkin(5) * t14;
t24 = t17 * t14;
t23 = -t19 * pkin(3) - t17 * qJ(4);
t22 = pkin(3) * t17 - qJ(4) * t19;
t21 = t23 * qJD(3) + t19 * qJD(4);
t13 = -0.2e1 * t24;
t12 = 0.2e1 * t24;
t9 = (-t15 + t16) * qJD(3);
t8 = -pkin(2) + t23;
t7 = -t19 * t30 - t20 * t31;
t6 = -t20 * t14 + t17 * t30;
t5 = t18 * t14 + t17 * t29;
t3 = t18 * t31 - t19 * t29;
t2 = t22 * qJD(3) - t17 * qJD(4);
t1 = 0.2e1 * (-0.1e1 + t32) * t18 * t29;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t29, 0, 0, 0, 0, 0, 0, 0, 0, t7, t6, t4, -pkin(2) * t30 + t33, 0, 0, 0, 0, 0, 0, t7, t4, -t6, -t20 * t2 + t8 * t30 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0.2e1 * t9, 0, t13, 0, 0, t17 * t28, t19 * t28, 0, 0, t12, 0, -0.2e1 * t9, 0, 0, t13, -0.2e1 * t2 * t19 + 0.2e1 * t8 * t31, 0, -0.2e1 * t8 * t14 - 0.2e1 * t2 * t17, 0.2e1 * t8 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t3, t21 * t18 - t22 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, -t31, 0, -t25, t26, 0, 0, 0, t14, 0, 0, t31, 0, -t25, t21, -t26, t21 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, qJ(4) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
