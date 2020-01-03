% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:50
% EndTime: 2019-12-31 18:14:51
% DurationCPUTime: 0.45s
% Computational Cost: add. (407->66), mult. (856->106), div. (0->0), fcn. (713->4), ass. (0->45)
t32 = sin(pkin(7));
t34 = cos(qJ(3));
t52 = cos(pkin(7));
t44 = t52 * t34;
t33 = sin(qJ(3));
t51 = t33 * qJD(3);
t17 = -qJD(3) * t44 + t32 * t51;
t20 = t32 * t34 + t52 * t33;
t14 = t20 * t17;
t18 = t20 * qJD(3);
t21 = -t32 * t33 + t44;
t53 = t21 * t18;
t62 = 0.2e1 * t53 + 0.2e1 * t14;
t61 = (t17 * t32 + t52 * t18) * pkin(3);
t56 = -pkin(1) - pkin(6);
t48 = qJ(4) - t56;
t22 = t48 * t33;
t43 = t48 * t34;
t12 = -t32 * t22 + t52 * t43;
t13 = -t52 * t22 - t32 * t43;
t16 = -qJD(3) * t43 - t33 * qJD(4);
t35 = -t34 * qJD(4) + t48 * t51;
t8 = t32 * t16 - t52 * t35;
t9 = t52 * t16 + t32 * t35;
t39 = t12 * t18 - t13 * t17 + t9 * t20 - t8 * t21;
t59 = -0.2e1 * t14;
t58 = 2 * qJD(2);
t57 = 2 * qJD(5);
t28 = t33 * pkin(3) + qJ(2);
t50 = t34 * qJD(3);
t23 = pkin(3) * t50 + qJD(2);
t49 = qJ(2) * qJD(3);
t47 = t33 * t50;
t46 = t12 * t8 + t13 * t9;
t45 = t56 * qJD(3);
t41 = t21 * t17 + t18 * t20;
t25 = t32 * pkin(3) + qJ(5);
t27 = -t52 * pkin(3) - pkin(4);
t37 = t20 * qJD(5) - t25 * t17 + t18 * t27;
t36 = 0.2e1 * t39;
t31 = qJ(2) * t58;
t11 = -0.2e1 * t53;
t10 = t20 * pkin(4) - t21 * qJ(5) + t28;
t5 = -t17 * pkin(4) + t18 * qJ(5) - t21 * qJD(5) + t23;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t31, -0.2e1 * t47, 0.2e1 * (t33 ^ 2 - t34 ^ 2) * qJD(3), 0, 0.2e1 * t47, 0, 0, 0.2e1 * qJD(2) * t33 + 0.2e1 * t34 * t49, 0.2e1 * qJD(2) * t34 - 0.2e1 * t33 * t49, 0, t31, t11, 0.2e1 * t41, 0, t59, 0, 0, -0.2e1 * t28 * t17 + 0.2e1 * t23 * t20, -0.2e1 * t28 * t18 + 0.2e1 * t23 * t21, -t36, 0.2e1 * t28 * t23 + 0.2e1 * t46, t11, 0, -0.2e1 * t41, 0, 0, t59, -0.2e1 * t10 * t17 + 0.2e1 * t5 * t20, -t36, 0.2e1 * t10 * t18 - 0.2e1 * t5 * t21, 0.2e1 * t10 * t5 + 0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t39, 0, 0, 0, 0, 0, 0, 0, t62, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, 0, -t50, 0, -t33 * t45, -t34 * t45, 0, 0, 0, 0, -t18, 0, t17, 0, -t8, -t9, t61, (t32 * t9 - t52 * t8) * pkin(3), 0, -t18, 0, 0, -t17, 0, -t8, -t37, t9, t13 * qJD(5) + t9 * t25 + t8 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0, -t61, 0, 0, 0, 0, 0, 0, -t18, 0, -t17, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t25 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, t23, 0, 0, 0, 0, 0, 0, -t17, 0, t18, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
