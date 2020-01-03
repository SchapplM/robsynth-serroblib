% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:12
% EndTime: 2019-12-31 17:41:13
% DurationCPUTime: 0.35s
% Computational Cost: add. (104->46), mult. (278->77), div. (0->0), fcn. (161->2), ass. (0->32)
t18 = cos(qJ(3));
t17 = sin(qJ(3));
t30 = t17 * qJ(4);
t33 = pkin(3) + pkin(4);
t34 = t33 * t18 + t30;
t6 = 0.2e1 * (-t17 ^ 2 + t18 ^ 2) * qJD(3);
t20 = 2 * qJD(4);
t32 = pkin(6) - qJ(5);
t31 = qJ(4) * t18;
t14 = t17 * qJD(3);
t29 = t17 * qJD(4);
t15 = t18 * qJD(3);
t28 = t18 * qJD(4);
t27 = -0.2e1 * pkin(2) * qJD(3);
t26 = pkin(6) * t14;
t25 = pkin(6) * t15;
t24 = t17 * t15;
t10 = t32 * t18;
t23 = -pkin(3) * t14 + t29;
t22 = -t18 * pkin(3) - t30;
t21 = t22 * qJD(3) + t28;
t16 = qJ(4) * t20;
t12 = -0.2e1 * t24;
t11 = 0.2e1 * t24;
t9 = t32 * t17;
t7 = -pkin(2) + t22;
t5 = pkin(2) + t34;
t4 = qJ(4) * t15 + t23;
t3 = qJD(3) * t10 - t17 * qJD(5);
t2 = -t18 * qJD(5) - t32 * t14;
t1 = (-pkin(4) * t17 + t31) * qJD(3) + t23;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t2 - t18 * t3 + (t10 * t18 + t17 * t9) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t6, 0, t12, 0, 0, t17 * t27, t18 * t27, 0, 0, t11, 0, -t6, 0, 0, t12, 0.2e1 * t7 * t14 + 0.2e1 * t4 * t18, 0, -0.2e1 * t7 * t15 + 0.2e1 * t4 * t17, -0.2e1 * t7 * t4, t11, -t6, 0, t12, 0, 0, 0.2e1 * t1 * t18 - 0.2e1 * t5 * t14, 0.2e1 * t1 * t17 + 0.2e1 * t5 * t15, -0.2e1 * t3 * t17 - 0.2e1 * t2 * t18 + 0.2e1 * (t10 * t17 - t18 * t9) * qJD(3), 0.2e1 * t5 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t9 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, t15, t4, 0, 0, 0, 0, 0, 0, -t14, t15, 0, t29 + (-t17 * t33 + t31) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t14, 0, -t25, t26, 0, 0, 0, t15, 0, 0, t14, 0, -t25, t21, -t26, t21 * pkin(6), 0, 0, -t15, 0, -t14, 0, -t3, t2, t34 * qJD(3) - t28, t2 * qJ(4) + t10 * qJD(4) - t3 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t16, 0, 0, 0, 0, 0, 0, 0, t20, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t15, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
