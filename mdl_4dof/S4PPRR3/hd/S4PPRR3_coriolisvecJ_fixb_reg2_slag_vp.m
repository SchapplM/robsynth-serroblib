% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PPRR3
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:27
% EndTime: 2019-12-31 16:17:28
% DurationCPUTime: 0.17s
% Computational Cost: add. (121->36), mult. (346->65), div. (0->0), fcn. (201->4), ass. (0->39)
t27 = (qJD(3) * qJD(4));
t40 = -2 * t27;
t13 = sin(qJ(3));
t16 = qJD(4) ^ 2;
t17 = qJD(3) ^ 2;
t39 = (t16 + t17) * t13;
t15 = cos(qJ(3));
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t5 = qJD(3) * pkin(5) + t13 * qJD(2);
t3 = -t14 * qJD(1) - t12 * t5;
t29 = t12 * qJD(1);
t37 = t14 * t5;
t4 = -t29 + t37;
t21 = t12 * t3 - t14 * t4;
t38 = t21 * t15;
t36 = t12 * t14;
t35 = t16 * t12;
t9 = t16 * t14;
t34 = t17 * t15;
t10 = t12 ^ 2;
t11 = t14 ^ 2;
t33 = t10 - t11;
t32 = t10 + t11;
t30 = qJD(3) * pkin(3);
t28 = t15 * qJD(2);
t26 = t17 * t36;
t25 = qJD(3) * t28;
t6 = -t28 - t30;
t24 = -t6 - t28;
t23 = t15 * t40;
t22 = t27 * t36;
t20 = qJD(3) * t24;
t19 = qJD(4) * (-t24 - t30);
t1 = t3 * qJD(4) + t14 * t25;
t7 = qJD(4) * t29;
t2 = -qJD(4) * t37 - t12 * t25 + t7;
t18 = t1 * t14 - t2 * t12 + (-t12 * t4 - t14 * t3) * qJD(4);
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t9, 0, t21 * qJD(4) - t1 * t12 - t2 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * t13, -t34, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t23 - t14 * t39, t12 * t39 + t14 * t23, t32 * t34, -qJD(3) * t38 + ((t6 - t28) * qJD(3) + t18) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t22, t33 * t40, t9, -0.2e1 * t22, -t35, 0, -pkin(5) * t9 + t12 * t19, pkin(5) * t35 + t14 * t19, -t32 * t25 + t18, (t38 + (-t6 - t30) * t13) * qJD(2) + t18 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t33 * t17, 0, t26, 0, 0, t7 + (t4 - t37) * qJD(4) + t12 * t20, t14 * t20, 0, 0;];
tauc_reg = t8;
