% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRR9
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
% MMD_reg [((4+1)*4/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:28
% EndTime: 2019-12-31 16:56:29
% DurationCPUTime: 0.22s
% Computational Cost: add. (103->49), mult. (303->110), div. (0->0), fcn. (197->4), ass. (0->50)
t50 = 2 * qJD(2);
t16 = cos(qJ(3));
t49 = t16 * pkin(6);
t15 = cos(qJ(4));
t11 = t15 ^ 2;
t13 = sin(qJ(4));
t48 = -t13 ^ 2 + t11;
t14 = sin(qJ(3));
t17 = -pkin(1) - pkin(5);
t47 = t14 * t17;
t46 = t16 * t17;
t10 = t14 ^ 2;
t12 = t16 ^ 2;
t45 = t10 - t12;
t44 = t10 + t12;
t43 = qJD(3) * t14;
t42 = qJD(3) * t15;
t41 = qJD(3) * t16;
t40 = qJD(3) * t17;
t39 = qJD(4) * t12;
t38 = qJD(4) * t13;
t37 = qJD(4) * t15;
t36 = qJD(4) * t16;
t35 = qJ(2) * qJD(3);
t34 = -0.2e1 * pkin(3) * qJD(4);
t33 = t15 * t47;
t32 = t13 * t36;
t31 = t15 * t36;
t30 = t13 * t41;
t29 = t13 * t37;
t28 = t14 * t42;
t27 = t15 * t41;
t26 = t14 * t41;
t25 = t48 * qJD(4);
t24 = t45 * qJD(3);
t23 = t17 * t30;
t22 = t13 * t28;
t21 = t17 * t27;
t20 = pkin(3) * t16 + pkin(6) * t14;
t19 = t14 * pkin(3) - t49;
t8 = qJ(2) + t19;
t18 = -t17 * t39 + t8 * t41;
t7 = t20 * qJD(3) + qJD(2);
t6 = t13 * t43 - t31;
t5 = t14 * t37 + t30;
t4 = -t28 - t32;
t3 = t14 * t38 - t27;
t2 = -t23 + t15 * t7 + (-t13 * t8 - t33) * qJD(4);
t1 = -t21 - t13 * t7 + (t13 * t47 - t15 * t8) * qJD(4);
t9 = [0, 0, 0, 0, t50, qJ(2) * t50, -0.2e1 * t26, 0.2e1 * t24, 0, 0, 0, 0.2e1 * qJD(2) * t14 + 0.2e1 * t16 * t35, 0.2e1 * qJD(2) * t16 - 0.2e1 * t14 * t35, -0.2e1 * t11 * t26 - 0.2e1 * t12 * t29, 0.4e1 * t16 * t22 - 0.2e1 * t48 * t39, -0.2e1 * t14 * t32 - 0.2e1 * t45 * t42, 0.2e1 * t13 * t24 - 0.2e1 * t14 * t31, 0.2e1 * t26, 0.2e1 * t18 * t15 + 0.2e1 * (t2 + t23) * t14, 0.2e1 * (t1 + t21) * t14 - 0.2e1 * t18 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t37, t44 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t41, 0, -t14 * t40, -t16 * t40, t16 * t25 - t22, -0.4e1 * t16 * t29 - t48 * t43, t5, -t3, 0, (-t13 * t46 - t15 * t20) * qJD(4) + (t13 * t19 - t33) * qJD(3), (t13 * t20 - t15 * t46) * qJD(4) + (-t15 * t49 + (pkin(3) * t15 + t13 * t17) * t14) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t41, 0, 0, 0, 0, 0, t4, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t29, 0.2e1 * t25, 0, 0, 0, t13 * t34, t15 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t6, t41, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t38, 0, -pkin(6) * t37, pkin(6) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
