% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR7
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
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:53
% EndTime: 2019-12-31 17:59:54
% DurationCPUTime: 0.28s
% Computational Cost: add. (137->54), mult. (350->109), div. (0->0), fcn. (244->6), ass. (0->47)
t17 = sin(qJ(4));
t12 = t17 ^ 2;
t19 = cos(qJ(4));
t14 = t19 ^ 2;
t26 = (t12 - t14) * qJD(4);
t49 = 2 * qJD(3);
t16 = sin(qJ(5));
t9 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(6);
t48 = t16 * t9;
t47 = t17 * t9;
t18 = cos(qJ(5));
t46 = t18 * t19;
t13 = t18 ^ 2;
t45 = t16 ^ 2 - t13;
t43 = t12 + t14;
t42 = qJD(5) * t14;
t41 = qJD(5) * t16;
t40 = qJD(5) * t18;
t39 = qJD(5) * t19;
t38 = t17 * qJD(4);
t37 = t19 * qJD(4);
t36 = t18 * t47;
t35 = -0.2e1 * pkin(4) * qJD(5);
t34 = t16 * t39;
t33 = t18 * t39;
t32 = t16 * t37;
t31 = t16 * t40;
t30 = t18 * t38;
t29 = t18 * t37;
t28 = t17 * t37;
t10 = sin(pkin(8)) * pkin(1) + qJ(3);
t27 = t45 * qJD(5);
t25 = t9 * t32;
t24 = t9 * t29;
t23 = t16 * t30;
t22 = pkin(4) * t19 + pkin(7) * t17;
t21 = t17 * pkin(4) - t19 * pkin(7);
t7 = t10 + t21;
t20 = t7 * t37 - t9 * t42;
t8 = t22 * qJD(4) + qJD(3);
t6 = t16 * t38 - t33;
t5 = t17 * t40 + t32;
t4 = t30 + t34;
t3 = t17 * t41 - t29;
t2 = -t25 + t18 * t8 + (-t16 * t7 - t36) * qJD(5);
t1 = -t24 - t16 * t8 + (t16 * t47 - t18 * t7) * qJD(5);
t11 = [0, 0, 0, 0, 0, t49, t10 * t49, -0.2e1 * t28, 0.2e1 * t26, 0, 0, 0, 0.2e1 * qJD(3) * t17 + 0.2e1 * t10 * t37, 0.2e1 * qJD(3) * t19 - 0.2e1 * t10 * t38, -0.2e1 * t13 * t28 - 0.2e1 * t14 * t31, 0.4e1 * t19 * t23 + 0.2e1 * t45 * t42, -0.2e1 * t17 * t34 - 0.2e1 * t18 * t26, 0.2e1 * t16 * t26 - 0.2e1 * t17 * t33, 0.2e1 * t28, 0.2e1 * t20 * t18 + 0.2e1 * (t2 + t25) * t17, 0.2e1 * (t1 + t24) * t17 - 0.2e1 * t20 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t40, t43 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t37, 0, -t9 * t38, -t9 * t37, -t19 * t27 - t23, -0.4e1 * t19 * t31 + t45 * t38, t5, -t3, 0, (-t22 * t18 - t19 * t48) * qJD(5) + (t21 * t16 - t36) * qJD(4), (t22 * t16 - t9 * t46) * qJD(5) + (-pkin(7) * t46 + (pkin(4) * t18 + t48) * t17) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t38, 0, 0, 0, 0, 0, t3, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t37, 0, 0, 0, 0, 0, -t4, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t31, -0.2e1 * t27, 0, 0, 0, t16 * t35, t18 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t6, t37, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t41, 0, -pkin(7) * t40, pkin(7) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
