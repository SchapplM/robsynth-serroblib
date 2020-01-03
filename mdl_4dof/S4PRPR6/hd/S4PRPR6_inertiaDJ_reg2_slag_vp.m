% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRPR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:41
% EndTime: 2019-12-31 16:24:42
% DurationCPUTime: 0.22s
% Computational Cost: add. (145->39), mult. (443->86), div. (0->0), fcn. (395->6), ass. (0->35)
t25 = cos(pkin(7));
t28 = cos(qJ(4));
t38 = t28 * t25;
t24 = sin(pkin(7));
t26 = sin(qJ(4));
t39 = t26 * t24;
t12 = -t38 + t39;
t27 = sin(qJ(2));
t8 = t12 * t27;
t35 = qJD(4) * t28;
t9 = qJD(4) * t39 - t25 * t35;
t41 = -0.2e1 * t9;
t13 = t28 * t24 + t26 * t25;
t10 = t13 * qJD(4);
t40 = 0.2e1 * t10;
t37 = pkin(5) + qJ(3);
t36 = t24 ^ 2 + t25 ^ 2;
t21 = t27 * qJD(2);
t29 = cos(qJ(2));
t34 = t29 * qJD(2);
t33 = t27 * t34;
t32 = t36 * t29;
t31 = t36 * qJD(3);
t30 = 0.2e1 * t31;
t15 = t37 * t24;
t16 = t37 * t25;
t6 = -t26 * t15 + t28 * t16;
t20 = -t25 * pkin(3) - pkin(2);
t7 = t13 * t27;
t5 = -t28 * t15 - t26 * t16;
t4 = qJD(4) * t8 - t13 * t34;
t3 = t27 * t10 + t12 * t34;
t2 = -t13 * qJD(3) - t6 * qJD(4);
t1 = t15 * t35 - qJD(3) * t38 + (qJD(3) * t24 + qJD(4) * t16) * t26;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t36) * t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t8 * t3 - 0.2e1 * t7 * t4 - 0.2e1 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t34, 0, 0, 0, 0, 0, 0, 0, 0, -t25 * t21, t24 * t21, qJD(2) * t32, t27 * t31 + (-pkin(2) * t27 + qJ(3) * t32) * qJD(2), 0, 0, 0, 0, 0, 0, -t29 * t10 + t12 * t21, t13 * t21 + t29 * t9, t8 * t10 + t3 * t12 - t4 * t13 - t7 * t9, t8 * t1 - t7 * t2 + t20 * t21 - t3 * t6 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, qJ(3) * t30, t13 * t41, -0.2e1 * t13 * t10 + 0.2e1 * t9 * t12, 0, t12 * t40, 0, 0, t20 * t40, t20 * t41, 0.2e1 * t1 * t12 - 0.2e1 * t6 * t10 - 0.2e1 * t2 * t13 + 0.2e1 * t5 * t9, -0.2e1 * t6 * t1 + 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, -t10, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
