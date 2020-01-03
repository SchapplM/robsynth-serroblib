% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:49
% EndTime: 2019-12-31 16:39:49
% DurationCPUTime: 0.23s
% Computational Cost: add. (282->50), mult. (568->87), div. (0->0), fcn. (279->4), ass. (0->50)
t18 = sin(pkin(6));
t23 = qJD(4) ^ 2;
t24 = qJD(1) ^ 2;
t58 = t18 * (t23 + t24);
t19 = cos(pkin(6));
t20 = sin(qJ(4));
t21 = cos(qJ(4));
t22 = -pkin(1) - pkin(2);
t12 = t22 * qJD(1) + qJD(2);
t44 = qJD(1) * qJ(2);
t8 = t18 * t12 + t19 * t44;
t6 = -qJD(1) * pkin(5) + t8;
t3 = t21 * qJD(3) - t20 * t6;
t4 = t20 * qJD(3) + t21 * t6;
t30 = t20 * t3 - t21 * t4;
t57 = t30 * t19;
t43 = qJD(1) * qJD(2);
t38 = t19 * t43;
t1 = t3 * qJD(4) + t21 * t38;
t2 = -t4 * qJD(4) - t20 * t38;
t25 = -(t20 * t4 + t21 * t3) * qJD(4) + t1 * t21 - t2 * t20;
t56 = t19 * t12;
t55 = t19 * t24;
t54 = t20 * t21;
t53 = t23 * t20;
t52 = t23 * t21;
t51 = t19 * qJ(2) + t18 * t22;
t16 = t20 ^ 2;
t17 = t21 ^ 2;
t50 = t16 - t17;
t49 = t16 + t17;
t46 = t18 * qJ(2);
t29 = t19 * t22 - t46;
t47 = qJD(1) * (pkin(3) - t29);
t45 = qJD(2) * t19;
t42 = qJD(1) * qJD(4);
t41 = t24 * t54;
t40 = 0.2e1 * t43;
t39 = 0.2e1 * t42;
t5 = -t56 + (pkin(3) + t46) * qJD(1);
t37 = -t5 - t45;
t36 = t19 * t39;
t35 = t18 * t40;
t34 = t42 * t54;
t32 = (-t18 * t44 + t56) * t18 - t8 * t19;
t28 = qJD(1) * (t5 - t45);
t10 = -pkin(5) + t51;
t27 = -t10 * t23 + t35;
t26 = qJD(4) * (t37 - t47);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, qJ(2) * t40, 0, 0, 0, 0, 0, 0, t35, 0.2e1 * t38, 0, ((-t18 * t29 + t19 * t51) * qJD(1) - t32) * qJD(2), 0.2e1 * t34, -t50 * t39, -t52, -0.2e1 * t34, t53, 0, t20 * t26 + t27 * t21, -t27 * t20 + t21 * t26, -t49 * t38 - t25, t25 * t10 + (-t57 + (t5 + t47) * t18) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t24 * qJ(2), 0, 0, 0, 0, 0, 0, -t18 * t24, -t55, 0, t32 * qJD(1), 0, 0, 0, 0, 0, 0, t20 * t36 - t21 * t58, t20 * t58 + t21 * t36, t49 * t55, qJD(1) * t57 + (t37 * qJD(1) + t25) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52, 0, -t30 * qJD(4) + t1 * t20 + t2 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t50 * t24, 0, t41, 0, 0, t20 * t28, t21 * t28, 0, 0;];
tauc_reg = t7;
