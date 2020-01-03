% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:05
% EndTime: 2019-12-31 18:26:05
% DurationCPUTime: 0.18s
% Computational Cost: add. (221->46), mult. (418->87), div. (0->0), fcn. (318->6), ass. (0->39)
t44 = 2 * qJD(2);
t33 = cos(qJ(3));
t34 = -pkin(1) - pkin(2);
t40 = qJD(3) * t33;
t31 = sin(qJ(3));
t41 = qJD(3) * t31;
t10 = qJ(2) * t41 - t33 * qJD(2) - t34 * t40;
t20 = t33 * qJ(2) + t31 * t34;
t11 = t31 * qJD(2) + t20 * qJD(3);
t28 = sin(pkin(8));
t29 = cos(pkin(8));
t1 = -t28 * t10 + t29 * t11;
t30 = sin(qJ(5));
t43 = t1 * t30;
t32 = cos(qJ(5));
t42 = t1 * t32;
t19 = -t31 * qJ(2) + t33 * t34 - pkin(3);
t8 = t28 * t19 + t29 * t20;
t39 = qJD(5) * t30;
t38 = qJD(5) * t32;
t26 = -t29 * pkin(3) - pkin(4);
t37 = 0.2e1 * qJD(5) * t26;
t36 = t30 * t38;
t7 = t29 * t19 - t28 * t20;
t5 = pkin(4) - t7;
t35 = qJD(5) * (-t26 + t5);
t17 = t28 * t33 + t29 * t31;
t16 = t28 * t31 - t29 * t33;
t25 = t28 * pkin(3) + pkin(7);
t22 = 0.2e1 * t36;
t21 = (-t30 ^ 2 + t32 ^ 2) * qJD(5);
t18 = 0.2e1 * t21;
t13 = t16 * qJD(3);
t12 = t17 * qJD(3);
t6 = -pkin(7) + t8;
t4 = t12 * t32 - t16 * t39;
t3 = t12 * t30 + t16 * t38;
t2 = -t29 * t10 - t28 * t11;
t9 = [0, 0, 0, 0, t44, qJ(2) * t44, 0, 0.2e1 * t11, -0.2e1 * t10, -0.2e1 * t7 * t1 + 0.2e1 * t8 * t2, t22, t18, 0, 0, 0, -0.2e1 * t5 * t39 + 0.2e1 * t42, -0.2e1 * t5 * t38 - 0.2e1 * t43; 0, 0, 0, 0, 0, 0, 0, t41, t40, t1 * t16 - t7 * t12 - t8 * t13 + t2 * t17, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t16 * t12 - 0.2e1 * t17 * t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t11, t10, (-t1 * t29 + t2 * t28) * pkin(3), -0.2e1 * t36, -0.2e1 * t21, 0, 0, 0, t30 * t35 - t42, t32 * t35 + t43; 0, 0, 0, 0, 0, 0, 0, -t41, -t40, (-t12 * t29 - t13 * t28) * pkin(3), 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t18, 0, 0, 0, t30 * t37, t32 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t39, 0, -t30 * t2 - t6 * t38, -t32 * t2 + t6 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t13 - t17 * t38, t32 * t13 + t17 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, 0, -t25 * t38, t25 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
