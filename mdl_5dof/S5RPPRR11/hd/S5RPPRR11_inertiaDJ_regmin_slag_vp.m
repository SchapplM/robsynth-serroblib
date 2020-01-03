% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:53
% EndTime: 2019-12-31 18:05:54
% DurationCPUTime: 0.30s
% Computational Cost: add. (125->51), mult. (343->113), div. (0->0), fcn. (227->4), ass. (0->47)
t17 = sin(qJ(4));
t11 = t17 ^ 2;
t19 = cos(qJ(4));
t13 = t19 ^ 2;
t29 = (t11 - t13) * qJD(4);
t18 = cos(qJ(5));
t12 = t18 ^ 2;
t16 = sin(qJ(5));
t48 = t16 ^ 2 - t12;
t30 = t48 * qJD(5);
t14 = -pkin(6) + qJ(2);
t41 = t17 * qJD(4);
t23 = t19 * qJD(2) - t14 * t41;
t26 = pkin(4) * t19 + pkin(7) * t17;
t49 = qJD(5) * t26 - t23;
t15 = pkin(1) + qJ(3);
t46 = t11 + t13;
t45 = qJD(5) * t16;
t44 = qJD(5) * t17;
t9 = qJD(5) * t18;
t43 = qJD(5) * t19;
t42 = t13 * qJD(2);
t40 = t19 * qJD(4);
t39 = qJ(2) * qJD(2);
t38 = -0.2e1 * pkin(4) * qJD(5);
t37 = t16 * t43;
t36 = t18 * t43;
t35 = t16 * t40;
t34 = t16 * t9;
t33 = t18 * t41;
t32 = t18 * t40;
t31 = t17 * t40;
t28 = t18 * t31;
t27 = -qJD(4) * t26 + t14 * t44 - qJD(3);
t25 = t17 * pkin(4) - t19 * pkin(7);
t24 = -t17 * qJD(2) - t14 * t40;
t8 = t25 + t15;
t22 = -qJD(5) * t8 + t24;
t21 = qJD(4) * t25 - t14 * t43;
t20 = 0.2e1 * qJD(2);
t6 = t16 * t41 - t36;
t5 = t17 * t9 + t35;
t4 = -t33 - t37;
t3 = t16 * t44 - t32;
t2 = t16 * t22 - t18 * t27;
t1 = t16 * t27 + t18 * t22;
t7 = [0, 0, 0, 0, t20, 0.2e1 * t39, t20, 0.2e1 * qJD(3), 0.2e1 * t15 * qJD(3) + 0.2e1 * t39, -0.2e1 * t31, 0.2e1 * t29, 0, 0, 0, 0.2e1 * qJD(3) * t17 + 0.2e1 * t15 * t40, 0.2e1 * qJD(3) * t19 - 0.2e1 * t15 * t41, -0.2e1 * t12 * t31 - 0.2e1 * t13 * t34, 0.2e1 * t13 * t30 + 0.4e1 * t16 * t28, -0.2e1 * t17 * t37 - 0.2e1 * t18 * t29, 0.2e1 * t16 * t29 - 0.2e1 * t17 * t36, 0.2e1 * t31, 0.2e1 * t8 * t32 - 0.2e1 * t16 * t42 + 0.2e1 * t2 * t17 + 0.2e1 * (-t13 * t9 + t16 * t31) * t14, -0.2e1 * t8 * t35 - 0.2e1 * t18 * t42 + 0.2e1 * t1 * t17 + 0.2e1 * (t13 * t45 + t28) * t14; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, 0, 0, 0, -t40, t41, 0, 0, 0, 0, 0, t3, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46 * t9, t46 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t40, 0, t23, t24, -t16 * t33 - t19 * t30, -0.4e1 * t19 * t34 + t41 * t48, t5, -t3, 0, t21 * t16 - t49 * t18, t49 * t16 + t21 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t40, 0, 0, 0, 0, 0, t4, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t34, -0.2e1 * t30, 0, 0, 0, t16 * t38, t18 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t6, t40, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t45, 0, -pkin(7) * t9, pkin(7) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
