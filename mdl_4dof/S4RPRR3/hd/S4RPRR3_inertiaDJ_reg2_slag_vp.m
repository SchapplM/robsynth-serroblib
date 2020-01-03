% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPRR3
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
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:18
% EndTime: 2019-12-31 16:49:19
% DurationCPUTime: 0.27s
% Computational Cost: add. (317->49), mult. (716->99), div. (0->0), fcn. (603->6), ass. (0->37)
t45 = qJD(3) + qJD(4);
t44 = 2 * qJD(3);
t24 = cos(qJ(3));
t40 = cos(qJ(4));
t30 = t40 * t24;
t22 = sin(qJ(4));
t23 = sin(qJ(3));
t39 = t22 * t23;
t13 = -t30 + t39;
t14 = t22 * t24 + t40 * t23;
t8 = t45 * t14;
t43 = t13 * t8;
t29 = t40 * qJD(4);
t7 = -qJD(3) * t30 - t24 * t29 + t39 * t45;
t42 = t14 * t7;
t18 = sin(pkin(7)) * pkin(1) + pkin(5);
t41 = pkin(6) + t18;
t38 = qJD(4) * t22;
t37 = t23 * qJD(3);
t36 = t24 * qJD(3);
t19 = -cos(pkin(7)) * pkin(1) - pkin(2);
t35 = t19 * t44;
t34 = pkin(3) * t37;
t33 = pkin(3) * t38;
t32 = t23 * t36;
t31 = t22 * t41;
t28 = t23 * t31;
t27 = pkin(3) * t29;
t26 = t41 * t40;
t25 = t23 * t26;
t12 = t41 * t24;
t6 = t40 * t12 - t28;
t15 = -t24 * pkin(3) + t19;
t5 = -t22 * t12 - t25;
t2 = -t6 * qJD(4) + (-t24 * t26 + t28) * qJD(3);
t1 = t12 * t38 + t25 * t45 + t31 * t36;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t32, (-t23 ^ 2 + t24 ^ 2) * t44, 0, -0.2e1 * t32, 0, 0, t23 * t35, t24 * t35, 0, 0, -0.2e1 * t42, 0.2e1 * t13 * t7 - 0.2e1 * t14 * t8, 0, 0.2e1 * t43, 0, 0, 0.2e1 * t13 * t34 + 0.2e1 * t15 * t8, 0.2e1 * t14 * t34 - 0.2e1 * t15 * t7, 0.2e1 * t1 * t13 - 0.2e1 * t2 * t14 + 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t6 * t1 + 0.2e1 * t15 * t34 + 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t14 - t2 * t13 - t5 * t8 - t6 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t42 + 0.2e1 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, -t37, 0, -t18 * t36, t18 * t37, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, (t40 * t7 - t22 * t8 + (-t40 * t13 + t14 * t22) * qJD(4)) * pkin(3), (t40 * t2 - t1 * t22 + (-t22 * t5 + t40 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t36, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, (-t40 * t8 - t22 * t7 + (t13 * t22 + t40 * t14) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t33, -0.2e1 * t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
