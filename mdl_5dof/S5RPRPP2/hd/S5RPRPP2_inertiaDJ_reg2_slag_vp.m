% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:12
% EndTime: 2019-12-31 18:11:13
% DurationCPUTime: 0.37s
% Computational Cost: add. (138->49), mult. (312->79), div. (0->0), fcn. (195->4), ass. (0->36)
t21 = cos(qJ(3));
t37 = pkin(3) + pkin(4);
t39 = t37 * t21;
t38 = 2 * qJD(3);
t20 = sin(qJ(3));
t9 = (-t20 ^ 2 + t21 ^ 2) * t38;
t23 = 2 * qJD(4);
t36 = t21 * pkin(3);
t35 = qJ(4) * t21;
t34 = t20 * qJ(4);
t13 = sin(pkin(7)) * pkin(1) + pkin(6);
t33 = qJ(5) - t13;
t16 = t20 * qJD(3);
t32 = t20 * qJD(4);
t17 = t21 * qJD(3);
t31 = t21 * qJD(4);
t14 = -cos(pkin(7)) * pkin(1) - pkin(2);
t30 = t14 * t38;
t29 = t20 * t17;
t28 = t13 * t16;
t27 = t13 * t17;
t8 = t33 * t21;
t26 = -pkin(3) * t16 + t32;
t25 = -t14 + t34;
t24 = t31 + (-t34 - t36) * qJD(3);
t18 = qJ(4) * t23;
t12 = -0.2e1 * t29;
t11 = 0.2e1 * t29;
t7 = t33 * t20;
t6 = -t25 - t36;
t5 = qJ(4) * t17 + t26;
t4 = t25 + t39;
t3 = -qJD(3) * t8 - t20 * qJD(5);
t2 = -t21 * qJD(5) + t33 * t16;
t1 = (-pkin(4) * t20 + t35) * qJD(3) + t26;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t9, 0, t12, 0, 0, t20 * t30, t21 * t30, 0, 0, t11, 0, -t9, 0, 0, t12, 0.2e1 * t6 * t16 + 0.2e1 * t5 * t21, 0, -0.2e1 * t6 * t17 + 0.2e1 * t5 * t20, -0.2e1 * t6 * t5, t11, -t9, 0, t12, 0, 0, 0.2e1 * t1 * t21 - 0.2e1 * t4 * t16, 0.2e1 * t1 * t20 + 0.2e1 * t4 * t17, -0.2e1 * t2 * t21 - 0.2e1 * t3 * t20 + 0.2e1 * (-t20 * t8 + t21 * t7) * qJD(3), 0.2e1 * t4 * t1 - 0.2e1 * t8 * t2 - 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * t20 - t3 * t21 + (-t20 * t7 - t21 * t8) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t16, 0, -t27, t28, 0, 0, 0, t17, 0, 0, t16, 0, -t27, t24, -t28, t24 * t13, 0, 0, -t17, 0, -t16, 0, -t3, t2, -t31 + (t34 + t39) * qJD(3), t2 * qJ(4) - t8 * qJD(4) - t3 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, t17, t5, 0, 0, 0, 0, 0, 0, -t16, t17, 0, t32 + (-t20 * t37 + t35) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t18, 0, 0, 0, 0, 0, 0, 0, t23, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t17, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
