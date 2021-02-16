% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:29
% EndTime: 2021-01-15 11:04:30
% DurationCPUTime: 0.21s
% Computational Cost: add. (171->59), mult. (441->102), div. (0->0), fcn. (264->4), ass. (0->45)
t34 = sin(qJ(3));
t29 = t34 * qJD(3);
t28 = pkin(3) * t29;
t35 = sin(qJ(2));
t46 = pkin(1) * qJD(2);
t42 = t35 * t46;
t13 = t28 + t42;
t36 = cos(qJ(3));
t26 = -t36 * pkin(3) - pkin(2);
t37 = cos(qJ(2));
t50 = t37 * pkin(1);
t18 = t26 - t50;
t31 = t36 * qJD(3);
t51 = t13 * t34 + t18 * t31;
t49 = -qJ(4) - pkin(6);
t25 = -pkin(2) - t50;
t48 = t25 * t31 + t34 * t42;
t33 = t34 ^ 2;
t47 = qJD(3) * t33 * pkin(3) + t26 * t31;
t24 = t35 * pkin(1) + pkin(6);
t45 = qJ(4) + t24;
t44 = pkin(2) * t29;
t43 = pkin(2) * t31;
t41 = t37 * t46;
t40 = pkin(3) * t31;
t8 = t18 * t29;
t16 = t26 * t29;
t39 = t34 * t31;
t38 = t25 * t29 - t36 * t42;
t32 = t36 * qJ(4);
t30 = t36 * qJD(4);
t23 = 0.2e1 * t39;
t21 = t36 * t41;
t20 = t36 * pkin(6) + t32;
t19 = t49 * t34;
t12 = 0.2e1 * (t36 ^ 2 - t33) * qJD(3);
t11 = t36 * t24 + t32;
t10 = t45 * t34;
t6 = -t34 * qJD(4) + t49 * t31;
t5 = -t49 * t29 - t30;
t4 = t5 * t36;
t3 = (-qJD(4) - t41) * t34 - t45 * t31;
t2 = t45 * t29 - t21 - t30;
t1 = t2 * t36;
t7 = [0, 0, 0, 0, -0.2e1 * t42, -0.2e1 * t41, t23, t12, 0, 0, 0, 0.2e1 * t38, 0.2e1 * t48, -0.2e1 * t13 * t36 + 0.2e1 * t8, 0.2e1 * t51, -0.2e1 * t3 * t34 - 0.2e1 * t1 + 0.2e1 * (t10 * t36 - t11 * t34) * qJD(3), -0.2e1 * t10 * t3 - 0.2e1 * t11 * t2 + 0.2e1 * t18 * t13; 0, 0, 0, 0, -t42, -t41, t23, t12, 0, 0, 0, t38 - t44, -t43 + t48, t16 + t8 + (-t13 - t28) * t36, t47 + t51, -t1 - t4 + (-t3 - t6) * t34 + ((t10 - t19) * t36 + (-t11 - t20) * t34) * qJD(3), pkin(3) * t8 - t10 * t6 - t11 * t5 + t13 * t26 + t3 * t19 - t2 * t20; 0, 0, 0, 0, 0, 0, t23, t12, 0, 0, 0, -0.2e1 * t44, -0.2e1 * t43, -0.2e1 * pkin(3) * t39 + 0.2e1 * t16, 0.2e1 * t47, -0.2e1 * t6 * t34 - 0.2e1 * t4 + 0.2e1 * (-t19 * t36 - t20 * t34) * qJD(3), 0.2e1 * pkin(3) * t16 + 0.2e1 * t19 * t6 - 0.2e1 * t20 * t5; 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, 0, -t24 * t31 - t34 * t41, t24 * t29 - t21, t3, t2, -t40, t3 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, 0, -pkin(6) * t31, pkin(6) * t29, t6, t5, -t40, t6 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t31, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t31, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
