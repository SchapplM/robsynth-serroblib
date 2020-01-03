% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:38
% EndTime: 2019-12-31 17:56:40
% DurationCPUTime: 0.36s
% Computational Cost: add. (255->44), mult. (458->84), div. (0->0), fcn. (313->6), ass. (0->38)
t24 = sin(qJ(5));
t22 = t24 ^ 2;
t26 = cos(qJ(5));
t23 = t26 ^ 2;
t41 = t22 + t23;
t40 = 2 * qJD(3);
t25 = sin(qJ(4));
t18 = sin(pkin(8)) * pkin(1) + qJ(3);
t29 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(3);
t37 = cos(qJ(4));
t7 = t37 * t18 + t25 * t29;
t3 = t25 * qJD(3) + t7 * qJD(4);
t39 = t3 * t24;
t38 = t3 * t26;
t36 = qJD(4) * t25;
t19 = t24 * qJD(5);
t20 = t26 * qJD(5);
t35 = -0.2e1 * pkin(4) * qJD(5);
t34 = t3 * t37;
t33 = t24 * t20;
t27 = t37 * t29;
t2 = -t37 * qJD(3) - qJD(4) * t27 + t18 * t36;
t1 = t41 * t2;
t6 = -t25 * t18 + t27;
t4 = pkin(4) - t6;
t32 = qJD(5) * (pkin(4) + t4);
t31 = qJD(4) * t37;
t30 = qJD(5) * t37;
t28 = t41 * t37;
t17 = -0.2e1 * t33;
t16 = 0.2e1 * t33;
t14 = (-t22 + t23) * qJD(5);
t11 = 0.2e1 * t14;
t10 = t24 * t30 + t26 * t36;
t9 = t24 * t36 - t26 * t30;
t8 = t28 * qJD(4);
t5 = -pkin(7) + t7;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t18 * t40, 0, 0, 0, 0, 0, 0, 0.2e1 * t3, -0.2e1 * t2, 0, -0.2e1 * t7 * t2 - 0.2e1 * t6 * t3, t16, t11, 0, t17, 0, 0, -0.2e1 * t4 * t19 + 0.2e1 * t38, -0.2e1 * t4 * t20 - 0.2e1 * t39, 0.2e1 * t1, -0.2e1 * t5 * t1 + 0.2e1 * t4 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t31, 0, -t34 - t2 * t25 + (-t25 * t6 + t37 * t7) * qJD(4), 0, 0, 0, 0, 0, 0, t10, -t9, -t8, -t34 - t25 * t1 + (t25 * t4 + t28 * t5) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-t37 + t28) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2, 0, 0, t17, -0.2e1 * t14, 0, t16, 0, 0, t24 * t32 - t38, t26 * t32 + t39, -t1, -t3 * pkin(4) - pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t31, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, t8, (-pkin(4) * t25 + t28 * pkin(7)) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t11, 0, t17, 0, 0, t24 * t35, t26 * t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, t19, 0, t24 * t2 - t5 * t20, t5 * t19 + t26 * t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25 * t20 - t24 * t31, t25 * t19 - t26 * t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t19, 0, -pkin(7) * t20, pkin(7) * t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
