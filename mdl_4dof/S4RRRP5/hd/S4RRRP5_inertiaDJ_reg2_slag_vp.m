% Calculate inertial parameters regressor of joint inertia matrix time derivative for
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
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:06
% EndTime: 2019-12-31 17:17:07
% DurationCPUTime: 0.34s
% Computational Cost: add. (416->73), mult. (1056->127), div. (0->0), fcn. (826->4), ass. (0->45)
t55 = qJD(2) + qJD(3);
t33 = 2 * qJD(4);
t54 = -pkin(6) - pkin(5);
t53 = cos(qJ(3));
t30 = sin(qJ(3));
t31 = sin(qJ(2));
t52 = t30 * t31;
t32 = cos(qJ(2));
t51 = t30 * t32;
t50 = qJD(3) * t30;
t49 = t31 * qJD(2);
t48 = t32 * qJD(2);
t18 = t53 * t31 + t51;
t10 = t55 * t18;
t41 = t53 * t32;
t17 = -t41 + t52;
t47 = 0.2e1 * t17 * t10;
t46 = -0.2e1 * pkin(1) * qJD(2);
t45 = pkin(2) * t49;
t44 = pkin(2) * t50;
t43 = t31 * t48;
t28 = -t32 * pkin(2) - pkin(1);
t20 = t54 * t32;
t38 = t54 * t53;
t36 = t31 * t38;
t11 = -t30 * t20 - t36;
t12 = -t53 * t20 + t54 * t52;
t35 = qJD(2) * t38;
t40 = t54 * qJD(2);
t4 = -qJD(3) * t36 - t20 * t50 - t31 * t35 - t40 * t51;
t39 = t53 * qJD(3);
t5 = -t20 * t39 - t32 * t35 + (qJD(3) * t54 + t40) * t52;
t42 = t11 * t5 - t12 * t4;
t9 = -qJD(2) * t41 - t32 * t39 + t55 * t52;
t37 = t18 * t10 - t9 * t17;
t34 = -0.2e1 * t12 * t10 - 0.2e1 * t11 * t9 + 0.2e1 * t4 * t17 + 0.2e1 * t5 * t18;
t29 = pkin(2) * t39;
t27 = -t53 * pkin(2) - pkin(3);
t25 = t30 * pkin(2) + qJ(4);
t24 = -0.2e1 * t44;
t21 = t29 + qJD(4);
t8 = t17 * pkin(3) - t18 * qJ(4) + t28;
t7 = -0.2e1 * t18 * t9;
t2 = t10 * pkin(3) + t9 * qJ(4) - t18 * qJD(4) + t45;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t43, 0.2e1 * (-t31 ^ 2 + t32 ^ 2) * qJD(2), 0, -0.2e1 * t43, 0, 0, t31 * t46, t32 * t46, 0, 0, t7, -0.2e1 * t37, 0, t47, 0, 0, 0.2e1 * t28 * t10 + 0.2e1 * t17 * t45, 0.2e1 * t18 * t45 - 0.2e1 * t28 * t9, t34, 0.2e1 * t28 * t45 + 0.2e1 * t42, t7, 0, 0.2e1 * t37, 0, 0, t47, 0.2e1 * t8 * t10 + 0.2e1 * t2 * t17, t34, -0.2e1 * t2 * t18 + 0.2e1 * t8 * t9, 0.2e1 * t8 * t2 + 0.2e1 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t49, 0, -pkin(5) * t48, pkin(5) * t49, 0, 0, 0, 0, -t9, 0, -t10, 0, -t5, t4, (t53 * t9 - t10 * t30 + (-t53 * t17 + t18 * t30) * qJD(3)) * pkin(2), (-t53 * t5 - t30 * t4 + (t11 * t30 + t53 * t12) * qJD(3)) * pkin(2), 0, -t9, 0, 0, t10, 0, -t5, -t25 * t10 - t21 * t17 + t18 * t44 - t27 * t9, -t4, t11 * t44 + t12 * t21 - t4 * t25 + t5 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -0.2e1 * t29, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0.2e1 * t21, 0.2e1 * t25 * t21 + 0.2e1 * t27 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, -t10, 0, -t5, t4, 0, 0, 0, -t9, 0, 0, t10, 0, -t5, pkin(3) * t9 - t10 * qJ(4) - t17 * qJD(4), -t4, -t5 * pkin(3) - t4 * qJ(4) + t12 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t29, 0, 0, 0, 0, 0, 0, 0, 0, -t44, 0, t33 + t29, -pkin(3) * t44 + t21 * qJ(4) + t25 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, qJ(4) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
