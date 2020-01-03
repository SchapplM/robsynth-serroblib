% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:14
% EndTime: 2019-12-31 18:01:14
% DurationCPUTime: 0.19s
% Computational Cost: add. (165->38), mult. (320->71), div. (0->0), fcn. (256->6), ass. (0->37)
t40 = 2 * qJD(2);
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t11 = t25 * t20 + t23 * t21;
t26 = -pkin(1) - pkin(2);
t28 = -t20 * qJ(2) + t21 * t26;
t13 = -pkin(3) + t28;
t14 = t21 * qJ(2) + t20 * t26;
t27 = t23 * t13 + t25 * t14;
t2 = t11 * qJD(2) + t27 * qJD(4);
t22 = sin(qJ(5));
t39 = t2 * t22;
t24 = cos(qJ(5));
t38 = t2 * t24;
t37 = t23 * t20;
t36 = qJD(4) * t25;
t35 = qJD(5) * t22;
t34 = qJD(5) * t24;
t33 = t20 * qJD(2);
t32 = t21 * qJD(2);
t31 = -0.2e1 * pkin(4) * qJD(5);
t30 = t22 * t34;
t5 = -t25 * t13 + t23 * t14 + pkin(4);
t29 = qJD(5) * (pkin(4) + t5);
t16 = 0.2e1 * t30;
t15 = (-t22 ^ 2 + t24 ^ 2) * qJD(5);
t12 = 0.2e1 * t15;
t10 = -t25 * t21 + t37;
t9 = t11 * qJD(4);
t8 = qJD(4) * t37 - t21 * t36;
t6 = -pkin(7) + t27;
t4 = -t10 * t35 + t9 * t24;
t3 = t10 * t34 + t9 * t22;
t1 = -t25 * t32 - t13 * t36 + (qJD(4) * t14 + t33) * t23;
t7 = [0, 0, 0, 0, t40, qJ(2) * t40, 0.2e1 * t33, 0.2e1 * t32, (t14 * t21 - t28 * t20) * t40, 0, 0.2e1 * t2, -0.2e1 * t1, t16, t12, 0, 0, 0, -0.2e1 * t5 * t35 + 0.2e1 * t38, -0.2e1 * t5 * t34 - 0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -0.2e1 * t30, -0.2e1 * t15, 0, 0, 0, t22 * t29 - t38, t24 * t29 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t12, 0, 0, 0, t22 * t31, t24 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t35, 0, t22 * t1 - t6 * t34, t24 * t1 + t6 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t34 + t22 * t8, t11 * t35 + t24 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t35, 0, -pkin(7) * t34, pkin(7) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
