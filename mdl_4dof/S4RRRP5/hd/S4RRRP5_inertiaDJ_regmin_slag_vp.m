% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRP5
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
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:05
% EndTime: 2019-12-31 17:17:06
% DurationCPUTime: 0.23s
% Computational Cost: add. (295->64), mult. (740->112), div. (0->0), fcn. (566->4), ass. (0->39)
t45 = qJD(2) + qJD(3);
t28 = 2 * qJD(4);
t44 = -pkin(6) - pkin(5);
t43 = cos(qJ(3));
t25 = sin(qJ(3));
t26 = sin(qJ(2));
t42 = t25 * t26;
t27 = cos(qJ(2));
t41 = t25 * t27;
t40 = qJD(3) * t25;
t39 = t26 * qJD(2);
t38 = t27 * qJD(2);
t37 = -0.2e1 * pkin(1) * qJD(2);
t36 = pkin(2) * t39;
t35 = pkin(2) * t40;
t23 = -t27 * pkin(2) - pkin(1);
t34 = t43 * t27;
t33 = t44 * qJD(2);
t32 = t43 * qJD(3);
t31 = t44 * t43;
t30 = t26 * t31;
t29 = qJD(2) * t31;
t14 = t43 * t26 + t41;
t24 = pkin(2) * t32;
t22 = -t43 * pkin(2) - pkin(3);
t20 = t25 * pkin(2) + qJ(4);
t19 = -0.2e1 * t35;
t16 = t24 + qJD(4);
t15 = t44 * t27;
t13 = -t34 + t42;
t8 = -t43 * t15 + t44 * t42;
t7 = -t25 * t15 - t30;
t6 = t45 * t14;
t5 = -qJD(2) * t34 - t27 * t32 + t45 * t42;
t4 = t13 * pkin(3) - t14 * qJ(4) + t23;
t3 = -t15 * t32 - t27 * t29 + (qJD(3) * t44 + t33) * t42;
t2 = -qJD(3) * t30 - t15 * t40 - t26 * t29 - t33 * t41;
t1 = t6 * pkin(3) + t5 * qJ(4) - t14 * qJD(4) + t36;
t9 = [0, 0, 0, 0.2e1 * t26 * t38, 0.2e1 * (-t26 ^ 2 + t27 ^ 2) * qJD(2), 0, 0, 0, t26 * t37, t27 * t37, -0.2e1 * t14 * t5, 0.2e1 * t5 * t13 - 0.2e1 * t14 * t6, 0, 0, 0, 0.2e1 * t13 * t36 + 0.2e1 * t23 * t6, 0.2e1 * t14 * t36 - 0.2e1 * t23 * t5, 0.2e1 * t1 * t13 + 0.2e1 * t4 * t6, 0.2e1 * t2 * t13 + 0.2e1 * t3 * t14 - 0.2e1 * t7 * t5 - 0.2e1 * t8 * t6, -0.2e1 * t1 * t14 + 0.2e1 * t4 * t5, 0.2e1 * t4 * t1 - 0.2e1 * t8 * t2 + 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, t38, -t39, 0, -pkin(5) * t38, pkin(5) * t39, 0, 0, -t5, -t6, 0, -t3, t2, -t3, -t16 * t13 + t14 * t35 - t20 * t6 - t22 * t5, -t2, t8 * t16 - t2 * t20 + t3 * t22 + t7 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -0.2e1 * t24, t19, 0, 0.2e1 * t16, 0.2e1 * t20 * t16 + 0.2e1 * t22 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, -t3, t2, -t3, pkin(3) * t5 - t6 * qJ(4) - t13 * qJD(4), -t2, -t3 * pkin(3) - t2 * qJ(4) + t8 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t24, -t35, 0, t28 + t24, -pkin(3) * t35 + t16 * qJ(4) + t20 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, qJ(4) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
