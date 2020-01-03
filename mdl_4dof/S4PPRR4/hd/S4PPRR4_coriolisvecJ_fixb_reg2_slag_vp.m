% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:40
% EndTime: 2019-12-31 16:18:40
% DurationCPUTime: 0.21s
% Computational Cost: add. (229->45), mult. (680->83), div. (0->0), fcn. (514->6), ass. (0->43)
t17 = sin(pkin(7));
t18 = cos(pkin(7));
t20 = sin(qJ(3));
t22 = cos(qJ(3));
t13 = t20 * t17 - t22 * t18;
t44 = t13 * qJD(1);
t14 = t22 * t17 + t20 * t18;
t10 = t14 * qJD(1);
t12 = t14 * qJD(3);
t8 = qJD(1) * t12;
t43 = t8 * t13;
t19 = sin(qJ(4));
t21 = cos(qJ(4));
t42 = t19 * t21;
t23 = qJD(4) ^ 2;
t41 = t23 * t19;
t40 = t23 * t21;
t15 = t19 ^ 2;
t16 = t21 ^ 2;
t39 = t15 - t16;
t38 = t15 + t16;
t37 = qJD(3) * pkin(3);
t36 = qJD(4) * t19;
t35 = qJD(4) * t21;
t11 = t13 * qJD(3);
t34 = t11 * qJD(3);
t33 = qJD(3) * qJD(4);
t24 = qJD(3) ^ 2;
t32 = t24 * t42;
t5 = t44 - t37;
t7 = qJD(1) * t11;
t31 = -qJD(3) * t5 + t7;
t30 = t33 * t42;
t6 = qJD(3) * pkin(5) + t10;
t3 = t21 * qJD(2) - t19 * t6;
t4 = t19 * qJD(2) + t21 * t6;
t29 = t19 * t3 - t21 * t4;
t28 = pkin(5) * t23 - t10 * qJD(3) + t8;
t27 = qJD(4) * (t5 - t44 - t37);
t1 = t3 * qJD(4) - t21 * t7;
t2 = -t4 * qJD(4) + t19 * t7;
t25 = t1 * t21 - t2 * t19 + (-t19 * t4 - t21 * t3) * qJD(4);
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * qJD(3), t34, 0, -t10 * t11 + t12 * t44 - t7 * t14 + t43, 0, 0, 0, 0, 0, 0, t11 * t36 - t14 * t40 + (-t12 * t21 + t13 * t36) * qJD(3), t11 * t35 + t14 * t41 + (t12 * t19 + t13 * t35) * qJD(3), -t38 * t34, t29 * t11 + t5 * t12 + t25 * t14 + t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t40, 0, -t29 * qJD(4) + t1 * t19 + t2 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t30, -0.2e1 * t39 * t33, t40, -0.2e1 * t30, -t41, 0, t19 * t27 - t28 * t21, t28 * t19 + t21 * t27, qJD(3) * t38 * t44 + t25, -t8 * pkin(3) + t25 * pkin(5) - t5 * t10 - t29 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t39 * t24, 0, t32, 0, 0, t31 * t19, t31 * t21, 0, 0;];
tauc_reg = t9;
