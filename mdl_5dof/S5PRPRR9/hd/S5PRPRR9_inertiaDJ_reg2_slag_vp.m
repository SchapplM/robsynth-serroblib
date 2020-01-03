% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:45
% EndTime: 2019-12-31 17:39:46
% DurationCPUTime: 0.36s
% Computational Cost: add. (168->42), mult. (371->82), div. (0->0), fcn. (226->4), ass. (0->37)
t23 = sin(qJ(5));
t21 = t23 ^ 2;
t25 = cos(qJ(5));
t22 = t25 ^ 2;
t40 = t21 + t22;
t39 = 2 * qJD(3);
t38 = -pkin(2) - pkin(3);
t24 = sin(qJ(4));
t35 = cos(qJ(4));
t11 = t35 * qJ(3) + t24 * t38;
t3 = t24 * qJD(3) + t11 * qJD(4);
t37 = t3 * t23;
t36 = t3 * t25;
t34 = qJD(4) * t24;
t17 = t23 * qJD(5);
t18 = t25 * qJD(5);
t33 = -0.2e1 * pkin(4) * qJD(5);
t32 = t3 * t35;
t31 = t23 * t18;
t27 = t35 * t38;
t2 = qJ(3) * t34 - t35 * qJD(3) - qJD(4) * t27;
t1 = t40 * t2;
t10 = -t24 * qJ(3) + t27;
t8 = pkin(4) - t10;
t30 = qJD(5) * (pkin(4) + t8);
t29 = qJD(4) * t35;
t28 = qJD(5) * t35;
t26 = t40 * t35;
t14 = -0.2e1 * t31;
t13 = 0.2e1 * t31;
t12 = (-t21 + t22) * qJD(5);
t9 = -pkin(7) + t11;
t7 = 0.2e1 * t12;
t6 = t23 * t28 + t25 * t34;
t5 = t23 * t34 - t25 * t28;
t4 = t26 * qJD(4);
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, qJ(3) * t39, 0, 0, 0, 0, 0, 0, 0.2e1 * t3, -0.2e1 * t2, 0, -0.2e1 * t10 * t3 - 0.2e1 * t11 * t2, t13, t7, 0, t14, 0, 0, -0.2e1 * t8 * t17 + 0.2e1 * t36, -0.2e1 * t8 * t18 - 0.2e1 * t37, 0.2e1 * t1, -0.2e1 * t9 * t1 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t29, 0, -t32 - t2 * t24 + (-t10 * t24 + t35 * t11) * qJD(4), 0, 0, 0, 0, 0, 0, t6, -t5, -t4, -t32 - t24 * t1 + (t24 * t8 + t26 * t9) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-t35 + t26) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2, 0, 0, t14, -0.2e1 * t12, 0, t13, 0, 0, t23 * t30 - t36, t25 * t30 + t37, -t1, -t3 * pkin(4) - pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t29, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, t4, (-pkin(4) * t24 + t26 * pkin(7)) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t7, 0, t14, 0, 0, t23 * t33, t25 * t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, t17, 0, -t9 * t18 + t23 * t2, t9 * t17 + t25 * t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t18 - t23 * t29, t24 * t17 - t25 * t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t17, 0, -pkin(7) * t18, pkin(7) * t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
