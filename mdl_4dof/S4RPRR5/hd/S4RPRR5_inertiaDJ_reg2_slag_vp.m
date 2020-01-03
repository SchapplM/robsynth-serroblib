% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:40
% EndTime: 2019-12-31 16:51:41
% DurationCPUTime: 0.26s
% Computational Cost: add. (168->42), mult. (369->82), div. (0->0), fcn. (224->4), ass. (0->37)
t21 = sin(qJ(4));
t19 = t21 ^ 2;
t23 = cos(qJ(4));
t20 = t23 ^ 2;
t40 = t19 + t20;
t39 = 2 * qJD(2);
t38 = -pkin(1) - pkin(2);
t22 = sin(qJ(3));
t35 = cos(qJ(3));
t11 = t35 * qJ(2) + t22 * t38;
t3 = t22 * qJD(2) + t11 * qJD(3);
t37 = t3 * t21;
t36 = t3 * t23;
t34 = qJD(3) * t22;
t33 = t21 * qJD(4);
t32 = t23 * qJD(4);
t31 = -0.2e1 * pkin(3) * qJD(4);
t30 = t3 * t35;
t29 = t21 * t32;
t25 = t35 * t38;
t2 = qJ(2) * t34 - t35 * qJD(2) - qJD(3) * t25;
t1 = t40 * t2;
t10 = -t22 * qJ(2) + t25;
t8 = pkin(3) - t10;
t28 = qJD(4) * (pkin(3) + t8);
t27 = qJD(3) * t35;
t26 = qJD(4) * t35;
t24 = t40 * t35;
t14 = -0.2e1 * t29;
t13 = 0.2e1 * t29;
t12 = (-t19 + t20) * qJD(4);
t9 = -pkin(6) + t11;
t7 = 0.2e1 * t12;
t6 = t21 * t26 + t23 * t34;
t5 = t21 * t34 - t23 * t26;
t4 = t24 * qJD(3);
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, qJ(2) * t39, 0, 0, 0, 0, 0, 0, 0.2e1 * t3, -0.2e1 * t2, 0, -0.2e1 * t10 * t3 - 0.2e1 * t11 * t2, t13, t7, 0, t14, 0, 0, -0.2e1 * t8 * t33 + 0.2e1 * t36, -0.2e1 * t8 * t32 - 0.2e1 * t37, 0.2e1 * t1, -0.2e1 * t9 * t1 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t27, 0, -t30 - t2 * t22 + (-t10 * t22 + t35 * t11) * qJD(3), 0, 0, 0, 0, 0, 0, t6, -t5, -t4, -t30 - t22 * t1 + (t22 * t8 + t24 * t9) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-t35 + t24) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2, 0, 0, t14, -0.2e1 * t12, 0, t13, 0, 0, t21 * t28 - t36, t23 * t28 + t37, -t1, -t3 * pkin(3) - pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t27, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, t4, (-pkin(3) * t22 + t24 * pkin(6)) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t7, 0, t14, 0, 0, t21 * t31, t23 * t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, t33, 0, t21 * t2 - t9 * t32, t23 * t2 + t9 * t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t27 - t22 * t32, t22 * t33 - t23 * t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t33, 0, -pkin(6) * t32, pkin(6) * t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
