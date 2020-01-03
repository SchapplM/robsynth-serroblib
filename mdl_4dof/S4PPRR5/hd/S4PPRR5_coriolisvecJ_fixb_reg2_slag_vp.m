% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:52
% EndTime: 2019-12-31 16:19:52
% DurationCPUTime: 0.17s
% Computational Cost: add. (122->39), mult. (363->75), div. (0->0), fcn. (208->4), ass. (0->39)
t9 = sin(qJ(4));
t7 = t9 ^ 2;
t11 = cos(qJ(4));
t8 = t11 ^ 2;
t39 = t7 - t8;
t38 = t7 + t8;
t37 = t11 * t9;
t10 = sin(qJ(3));
t12 = cos(qJ(3));
t6 = t12 * qJD(1) + t10 * qJD(2);
t4 = t6 * qJD(3);
t36 = t4 * t10;
t35 = t4 * t12;
t14 = qJD(3) ^ 2;
t34 = t10 * t14;
t33 = t12 * t14;
t13 = qJD(4) ^ 2;
t32 = t13 + t14;
t31 = qJD(3) * pkin(3);
t30 = qJD(3) * qJD(4);
t29 = t14 * t37;
t17 = t10 * qJD(1) - t12 * qJD(2);
t28 = t38 * t17;
t3 = t17 * qJD(3);
t27 = t38 * t3;
t26 = t10 * t38;
t25 = t12 * t38;
t24 = 0.2e1 * t30;
t23 = t32 * t9;
t22 = t32 * t11;
t1 = t17 - t31;
t21 = -qJD(3) * t1 + t3;
t20 = t10 * t24;
t19 = -0.2e1 * t12 * t30;
t18 = t30 * t37;
t16 = pkin(5) * t13;
t15 = qJD(4) * (t1 - t17 - t31);
t2 = qJD(3) * pkin(5) + t6;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t34, 0, t36 - t3 * t12 + (-t10 * t6 + t12 * t17) * qJD(3), 0, 0, 0, 0, 0, 0, -t12 * t22 + t9 * t20, t11 * t20 + t12 * t23, -t14 * t26, t36 - t3 * t25 + (t1 * t12 - t2 * t26) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t33, 0, -t3 * t10 - t35 + (t10 * t17 + t12 * t6) * qJD(3), 0, 0, 0, 0, 0, 0, -t10 * t22 + t9 * t19, t10 * t23 + t11 * t19, t14 * t25, -t35 - t10 * t27 + (t1 * t10 + t2 * t25) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t18, -t39 * t24, t13 * t11, -0.2e1 * t18, -t13 * t9, 0, -t16 * t11 + t9 * t15, t11 * t15 + t16 * t9, qJD(3) * t28 - t27, -t4 * pkin(3) - pkin(5) * t27 - t1 * t6 + t2 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t39 * t14, 0, t29, 0, 0, t21 * t9, t21 * t11, 0, 0;];
tauc_reg = t5;
