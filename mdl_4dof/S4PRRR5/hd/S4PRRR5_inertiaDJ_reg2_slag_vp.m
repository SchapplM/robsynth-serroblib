% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:46
% EndTime: 2019-12-31 16:33:47
% DurationCPUTime: 0.28s
% Computational Cost: add. (154->39), mult. (484->78), div. (0->0), fcn. (393->6), ass. (0->38)
t24 = sin(qJ(3));
t25 = sin(qJ(2));
t40 = cos(qJ(3));
t41 = cos(qJ(2));
t8 = t24 * t25 - t40 * t41;
t23 = sin(qJ(4));
t21 = t23 ^ 2;
t26 = cos(qJ(4));
t22 = t26 ^ 2;
t45 = t21 + t22;
t44 = qJD(2) + qJD(3);
t9 = t24 * t41 + t40 * t25;
t5 = t44 * t9;
t43 = t8 * t5;
t42 = t24 * t8;
t32 = t40 * pkin(2);
t19 = -t32 - pkin(3);
t20 = t26 * qJD(4);
t37 = pkin(2) * qJD(3);
t33 = t24 * t37;
t38 = t19 * t20 + t23 * t33;
t36 = t23 * qJD(4);
t35 = pkin(3) * t36;
t34 = pkin(3) * t20;
t31 = t23 * t20;
t4 = t8 * t44;
t1 = t45 * t4;
t30 = qJD(3) * t32;
t28 = t19 * t36 - t26 * t33;
t27 = t45 * t40;
t18 = t24 * pkin(2) + pkin(6);
t14 = -0.2e1 * t31;
t13 = 0.2e1 * t31;
t7 = 0.2e1 * (-t21 + t22) * qJD(4);
t6 = t27 * t37;
t3 = -t5 * t26 + t8 * t36;
t2 = t8 * t20 + t5 * t23;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t9 * t4 + 0.2e1 * t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t9 * t1 + 0.2e1 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25 * qJD(2), -t41 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, 0, (-t40 * t5 - t24 * t4 + (t40 * t9 + t42) * qJD(3)) * pkin(2), 0, 0, 0, 0, 0, 0, t3, t2, -t1, t5 * t19 - t18 * t1 + (t27 * t9 + t42) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t33, -0.2e1 * t30, 0, 0, t13, t7, 0, t14, 0, 0, 0.2e1 * t28, 0.2e1 * t38, 0.2e1 * t6, 0.2e1 * (t27 * t18 + t19 * t24) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3, t2, -t1, -t5 * pkin(3) - pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t30, 0, 0, t13, t7, 0, t14, 0, 0, t28 - t35, -t34 + t38, t6, (-pkin(3) * t24 + t27 * pkin(6)) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t7, 0, t14, 0, 0, -0.2e1 * t35, -0.2e1 * t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t20 + t23 * t4, t26 * t4 + t9 * t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t36, 0, -t18 * t20 - t23 * t30, t18 * t36 - t26 * t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t36, 0, -pkin(6) * t20, pkin(6) * t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
