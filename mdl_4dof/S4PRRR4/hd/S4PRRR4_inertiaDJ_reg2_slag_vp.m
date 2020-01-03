% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4PRRR4
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:39
% EndTime: 2019-12-31 16:32:40
% DurationCPUTime: 0.26s
% Computational Cost: add. (246->47), mult. (645->97), div. (0->0), fcn. (532->4), ass. (0->34)
t41 = qJD(3) + qJD(4);
t40 = -pkin(6) - pkin(5);
t21 = cos(qJ(3));
t37 = cos(qJ(4));
t27 = t37 * t21;
t19 = sin(qJ(4));
t20 = sin(qJ(3));
t36 = t19 * t20;
t12 = -t27 + t36;
t13 = t19 * t21 + t37 * t20;
t6 = t41 * t13;
t39 = t12 * t6;
t26 = t37 * qJD(4);
t5 = -qJD(3) * t27 - t21 * t26 + t41 * t36;
t38 = t13 * t5;
t35 = qJD(4) * t19;
t34 = t20 * qJD(3);
t33 = t21 * qJD(3);
t32 = -0.2e1 * pkin(2) * qJD(3);
t31 = pkin(3) * t34;
t30 = pkin(3) * t35;
t29 = t19 * t40;
t28 = t20 * t33;
t25 = t20 * t29;
t24 = pkin(3) * t26;
t23 = t40 * t37;
t22 = t20 * t23;
t14 = t40 * t21;
t8 = -t37 * t14 + t25;
t18 = -t21 * pkin(3) - pkin(2);
t7 = t19 * t14 + t22;
t2 = -t8 * qJD(4) + (t21 * t23 - t25) * qJD(3);
t1 = -t14 * t35 - t41 * t22 - t29 * t33;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t38 + 0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t1 - t12 * t2 - t5 * t8 - t6 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t28, 0.2e1 * (-t20 ^ 2 + t21 ^ 2) * qJD(3), 0, -0.2e1 * t28, 0, 0, t20 * t32, t21 * t32, 0, 0, -0.2e1 * t38, 0.2e1 * t12 * t5 - 0.2e1 * t13 * t6, 0, 0.2e1 * t39, 0, 0, 0.2e1 * t12 * t31 + 0.2e1 * t18 * t6, 0.2e1 * t13 * t31 - 0.2e1 * t18 * t5, 0.2e1 * t1 * t12 - 0.2e1 * t2 * t13 + 0.2e1 * t7 * t5 - 0.2e1 * t8 * t6, -0.2e1 * t8 * t1 + 0.2e1 * t18 * t31 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t33, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, (-t37 * t6 - t19 * t5 + (t12 * t19 + t37 * t13) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t34, 0, -pkin(5) * t33, pkin(5) * t34, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, (t37 * t5 - t19 * t6 + (-t37 * t12 + t13 * t19) * qJD(4)) * pkin(3), (t37 * t2 - t1 * t19 + (-t19 * t7 + t37 * t8) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t30, -0.2e1 * t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
