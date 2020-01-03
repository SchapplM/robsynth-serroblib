% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:30
% EndTime: 2019-12-31 17:33:31
% DurationCPUTime: 0.26s
% Computational Cost: add. (245->58), mult. (491->83), div. (0->0), fcn. (256->4), ass. (0->49)
t13 = sin(qJ(5));
t15 = cos(qJ(5));
t17 = -pkin(3) - pkin(6);
t16 = cos(qJ(3));
t39 = t16 * qJD(2);
t31 = qJD(4) - t39;
t24 = -t17 * qJD(3) - t31;
t3 = t13 * qJD(1) - t15 * t24;
t14 = sin(qJ(3));
t41 = t14 * qJD(2);
t33 = qJD(3) * t41;
t1 = t3 * qJD(5) + t13 * t33;
t23 = t13 * t24;
t40 = t15 * qJD(1);
t51 = qJD(5) * t40 + t15 * t33;
t2 = qJD(5) * t23 + t51;
t4 = t23 + t40;
t21 = -t1 * t13 - t2 * t15 + (t13 * t3 + t15 * t4) * qJD(5);
t38 = qJD(3) * qJ(4);
t8 = t38 + t41;
t42 = t8 * qJD(3);
t55 = t21 + t42;
t32 = -t8 + t41;
t54 = qJD(5) * (-t32 + t38);
t53 = t32 * qJD(3);
t52 = t16 * t8;
t18 = qJD(5) ^ 2;
t50 = t18 * t13;
t49 = t18 * t15;
t19 = qJD(3) ^ 2;
t48 = t19 * t14;
t47 = t19 * t16;
t11 = t13 ^ 2;
t12 = t15 ^ 2;
t46 = t11 - t12;
t45 = t11 + t12;
t44 = t18 + t19;
t43 = qJD(3) * pkin(3);
t37 = qJD(3) * qJD(5);
t36 = t15 * t19 * t13;
t35 = t44 * t16;
t34 = t15 * t37;
t30 = t13 * t34;
t28 = t13 * t4 - t15 * t3;
t6 = (qJD(4) + t39) * qJD(3);
t26 = t6 * qJ(4) + t8 * qJD(4);
t22 = t31 * qJD(3) - t17 * t18 + t6;
t7 = t31 - t43;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t50, 0, -t28 * qJD(5) - t1 * t15 + t2 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t47, t6 * t14 + (t52 + (t7 - t39) * t14) * qJD(3), 0, 0, 0, 0, 0, 0, t13 * t35 + 0.2e1 * t14 * t34, -0.2e1 * t14 * t13 * t37 + t15 * t35, -t45 * t48, (-t28 * qJD(3) + t6) * t14 + t55 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * qJD(4), (-t52 + (-t7 - t43) * t14) * qJD(2) + t26, -0.2e1 * t30, 0.2e1 * t46 * t37, -t50, 0.2e1 * t30, -t49, 0, t22 * t13 + t15 * t54, -t13 * t54 + t22 * t15, t45 * t33 + t21, (t28 * t14 - t52) * qJD(2) - t21 * t17 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t53, 0, 0, 0, 0, 0, 0, -t44 * t13, -t44 * t15, 0, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t46 * t19, 0, -t36, 0, 0, -t15 * t42 + (t23 - t4) * qJD(5) + t51, -t13 * t53, 0, 0;];
tauc_reg = t5;
