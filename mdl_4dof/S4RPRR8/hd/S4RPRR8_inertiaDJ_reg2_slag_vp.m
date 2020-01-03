% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPRR8
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:16
% DurationCPUTime: 0.28s
% Computational Cost: add. (319->51), mult. (654->89), div. (0->0), fcn. (534->4), ass. (0->40)
t21 = cos(qJ(3));
t41 = sin(qJ(4));
t32 = t41 * t21;
t19 = sin(qJ(3));
t20 = cos(qJ(4));
t40 = t20 * t19;
t10 = t32 + t40;
t33 = t41 * t19;
t39 = t20 * t21;
t11 = -t33 + t39;
t47 = qJD(3) + qJD(4);
t3 = t47 * t10;
t28 = qJD(3) * t33;
t31 = qJD(4) * t41;
t4 = -t19 * t31 + t47 * t39 - t28;
t48 = ((-t10 * t20 + t41 * t11) * qJD(4) + t20 * t3 - t41 * t4) * pkin(3);
t45 = 2 * qJD(2);
t44 = t10 * t4;
t43 = t11 * t3;
t22 = -pkin(1) - pkin(5);
t42 = pkin(6) - t22;
t38 = t19 * qJD(3);
t37 = t21 * qJD(3);
t36 = qJ(2) * qJD(3);
t35 = qJD(4) * t20 * pkin(3);
t34 = t19 * t37;
t30 = t42 * t39;
t29 = pkin(3) * t31;
t27 = -t43 + t44;
t26 = t42 * t32;
t12 = t42 * t19;
t6 = -t20 * t12 - t26;
t1 = -t12 * t31 - t42 * t28 + t47 * t30;
t2 = -qJD(4) * t6 + (t42 * t40 + t26) * qJD(3);
t5 = t41 * t12 - t30;
t23 = t1 * t10 - t2 * t11 + t5 * t3 - t6 * t4;
t18 = qJ(2) * t45;
t16 = t19 * pkin(3) + qJ(2);
t13 = pkin(3) * t37 + qJD(2);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t18, -0.2e1 * t34, 0.2e1 * (t19 ^ 2 - t21 ^ 2) * qJD(3), 0, 0.2e1 * t34, 0, 0, 0.2e1 * qJD(2) * t19 + 0.2e1 * t21 * t36, 0.2e1 * qJD(2) * t21 - 0.2e1 * t19 * t36, 0, t18, -0.2e1 * t43, 0.2e1 * t3 * t10 - 0.2e1 * t11 * t4, 0, 0.2e1 * t44, 0, 0, 0.2e1 * t13 * t10 + 0.2e1 * t16 * t4, 0.2e1 * t13 * t11 - 0.2e1 * t16 * t3, 0.2e1 * t23, -0.2e1 * t6 * t1 + 0.2e1 * t16 * t13 + 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t27, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, -t37, 0, -t22 * t38, -t22 * t37, 0, 0, 0, 0, -t3, 0, -t4, 0, t2, t1, t48, (-t41 * t1 + t2 * t20 + (t20 * t6 - t41 * t5) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t37, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t29, -0.2e1 * t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, -t4, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
