% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:35
% EndTime: 2019-12-31 16:50:36
% DurationCPUTime: 0.20s
% Computational Cost: add. (101->43), mult. (312->99), div. (0->0), fcn. (215->6), ass. (0->48)
t17 = sin(qJ(3));
t50 = -0.4e1 * t17;
t16 = sin(qJ(4));
t9 = sin(pkin(7)) * pkin(1) + pkin(5);
t49 = t16 * t9;
t18 = cos(qJ(4));
t19 = cos(qJ(3));
t48 = t18 * t19;
t13 = t18 ^ 2;
t47 = t16 ^ 2 - t13;
t12 = t17 ^ 2;
t46 = -t19 ^ 2 + t12;
t45 = qJD(3) * t17;
t44 = qJD(3) * t18;
t43 = qJD(3) * t19;
t42 = qJD(4) * t12;
t41 = qJD(4) * t16;
t40 = qJD(4) * t18;
t39 = qJD(4) * t19;
t38 = t19 * t49;
t37 = t9 * t48;
t36 = -0.2e1 * pkin(3) * qJD(4);
t10 = -cos(pkin(7)) * pkin(1) - pkin(2);
t35 = 0.2e1 * qJD(3) * t10;
t34 = t16 * t39;
t33 = t18 * t39;
t32 = t16 * t45;
t31 = t16 * t40;
t30 = t17 * t43;
t29 = t17 * t44;
t28 = t18 * t43;
t27 = t47 * qJD(4);
t26 = t46 * qJD(3);
t25 = t9 * t32;
t24 = t9 * t29;
t23 = t16 * t28;
t22 = -t19 * pkin(3) - t17 * pkin(6);
t21 = pkin(3) * t17 - pkin(6) * t19;
t7 = t10 + t22;
t20 = t9 * t42 + t7 * t45;
t8 = t21 * qJD(3);
t6 = t32 - t33;
t5 = -t16 * t43 - t17 * t40;
t4 = t29 + t34;
t3 = t17 * t41 - t28;
t2 = t25 + t18 * t8 + (-t16 * t7 - t37) * qJD(4);
t1 = t24 - t16 * t8 + (-t18 * t7 + t38) * qJD(4);
t11 = [0, 0, 0, 0, 0.2e1 * t30, -0.2e1 * t26, 0, 0, 0, t17 * t35, t19 * t35, -0.2e1 * t12 * t31 + 0.2e1 * t13 * t30, t23 * t50 + 0.2e1 * t47 * t42, 0.2e1 * t17 * t34 + 0.2e1 * t46 * t44, -0.2e1 * t16 * t26 + 0.2e1 * t17 * t33, -0.2e1 * t30, 0.2e1 * (-t2 + t25) * t19 + 0.2e1 * t20 * t18, 0.2e1 * (-t1 + t24) * t19 - 0.2e1 * t20 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t43, -t45, 0, -t9 * t43, t9 * t45, -t17 * t27 + t23, t31 * t50 - t47 * t43, t6, t4, 0, (pkin(6) * t48 + (-pkin(3) * t18 + t49) * t17) * qJD(4) + (t22 * t16 - t37) * qJD(3), (t17 * t18 * t9 + t21 * t16) * qJD(4) + (t22 * t18 + t38) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t43, 0, 0, 0, 0, 0, -t4, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t31, -0.2e1 * t27, 0, 0, 0, t16 * t36, t18 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t5, t45, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t41, 0, -pkin(6) * t40, pkin(6) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
