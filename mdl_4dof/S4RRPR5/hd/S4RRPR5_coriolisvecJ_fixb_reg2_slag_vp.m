% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:33
% EndTime: 2019-12-31 17:03:34
% DurationCPUTime: 0.25s
% Computational Cost: add. (207->64), mult. (432->82), div. (0->0), fcn. (164->4), ass. (0->57)
t21 = sin(qJ(2));
t49 = pkin(1) * qJD(2);
t41 = t21 * t49;
t17 = qJD(1) + qJD(2);
t44 = qJD(1) + t17;
t4 = t44 * t41;
t43 = qJD(2) - t17;
t16 = t17 ^ 2;
t24 = -pkin(2) - pkin(6);
t50 = pkin(1) * qJD(1);
t39 = t21 * t50;
t48 = t17 * qJ(3);
t10 = t39 + t48;
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t46 = qJD(4) * t22;
t23 = cos(qJ(2));
t36 = qJD(1) * t49;
t33 = t23 * t36;
t45 = t17 * qJD(3);
t8 = t33 + t45;
t58 = t10 * t46 + t8 * t20;
t57 = t10 * t17;
t56 = t10 * t23;
t25 = qJD(4) ^ 2;
t55 = t25 * t20;
t54 = t25 * t22;
t53 = -t16 - t25;
t18 = t20 ^ 2;
t19 = t22 ^ 2;
t52 = t18 - t19;
t51 = t18 + t19;
t47 = qJD(4) * t20;
t42 = t22 * t16 * t20;
t40 = t23 * t49;
t38 = t23 * t50;
t37 = -t23 * pkin(1) - pkin(2);
t35 = t17 * t20 * t46;
t34 = t21 * t36;
t32 = t43 * t50;
t31 = qJD(3) - t38;
t13 = qJD(3) + t40;
t15 = t21 * pkin(1) + qJ(3);
t30 = t10 * t13 + t8 * t15;
t14 = -pkin(6) + t37;
t29 = t13 * t17 - t14 * t25;
t28 = t8 * qJ(3) + t10 * qJD(3);
t27 = t15 * t17 + t41;
t26 = t34 - t57;
t12 = -0.2e1 * t35;
t11 = 0.2e1 * t35;
t9 = -t17 * pkin(2) + t31;
t6 = t8 * t22;
t3 = t21 * t32;
t2 = t24 * t17 + t31;
t1 = 0.2e1 * t52 * t17 * qJD(4);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t44 * t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t33 + (qJD(3) + t13) * t17, (t37 * qJD(1) + t9) * t41 + t30, t12, t1, -t55, t11, -t54, 0, t29 * t20 + t27 * t46 + t58, t6 + t29 * t22 + (-t10 - t27) * t47, -t51 * t4, t30 + (qJD(1) * t14 + t2) * t41 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t43 * t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t23 * t32 + 0.2e1 * t45, (-t56 + (-pkin(2) * qJD(2) - t9) * t21) * t50 + t28, t12, t1, -t55, t11, -t54, 0, -t39 * t46 - t24 * t55 + (qJ(3) * t46 + t31 * t20) * t17 + t58, t6 + (t31 * t17 - t24 * t25) * t22 + (-t10 + t39 - t48) * t47, -t43 * t51 * t39, (-t56 + (qJD(2) * t24 - t2) * t21 * t51) * t50 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t26, 0, 0, 0, 0, 0, 0, t53 * t20, t53 * t22, 0, t51 * t34 - t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t52 * t16, 0, -t42, 0, 0, t26 * t22, -t26 * t20, 0, 0;];
tauc_reg = t5;
