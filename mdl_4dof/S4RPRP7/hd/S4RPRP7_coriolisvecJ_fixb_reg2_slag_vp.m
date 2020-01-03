% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:17
% EndTime: 2019-12-31 16:47:17
% DurationCPUTime: 0.23s
% Computational Cost: add. (174->53), mult. (398->82), div. (0->0), fcn. (148->2), ass. (0->47)
t20 = sin(qJ(3));
t22 = -pkin(1) - pkin(5);
t15 = t22 * qJD(1) + qJD(2);
t21 = cos(qJ(3));
t44 = t21 * t15;
t3 = (qJD(4) + t44) * qJD(3);
t39 = qJD(3) * pkin(3);
t30 = -qJD(4) + t39;
t4 = -t30 - t44;
t35 = qJD(3) * qJ(4);
t46 = t20 * t15;
t6 = t35 + t46;
t48 = t21 * t6;
t51 = ((-t4 + t44) * t20 - t48) * qJD(3) - t3 * t20;
t33 = 2 * qJD(1);
t49 = 0.2e1 * qJD(3);
t45 = t20 * t21;
t23 = qJD(3) ^ 2;
t43 = t23 * t20;
t42 = t23 * t21;
t19 = t21 ^ 2;
t41 = t20 ^ 2 - t19;
t24 = qJD(1) ^ 2;
t40 = t23 + t24;
t28 = pkin(3) * t21 + qJ(4) * t20;
t2 = t28 * qJD(3) - t21 * qJD(4) + qJD(2);
t1 = qJD(1) * t2;
t9 = t20 * pkin(3) - t21 * qJ(4) + qJ(2);
t5 = qJD(1) * t9;
t38 = t24 * qJ(2);
t37 = t5 * qJD(1);
t36 = qJ(2) * qJD(3);
t34 = qJD(1) * qJD(3);
t32 = qJD(2) * t33;
t29 = t34 * t45;
t27 = t5 * t49;
t26 = -t22 * t23 + 0.2e1 * t1;
t17 = qJ(2) * t32;
t16 = t24 * t45;
t14 = -0.2e1 * t29;
t13 = 0.2e1 * t29;
t12 = t40 * t21;
t11 = t40 * t20;
t10 = t41 * t24;
t8 = t28 * qJD(1);
t7 = t41 * t34;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t17, t14, 0.2e1 * t7, -t43, t13, -t42, 0, -t22 * t43 + (qJD(2) * t20 + t21 * t36) * t33, -t22 * t42 + (qJD(2) * t21 - t20 * t36) * t33, 0, t17, t14, -t43, -0.2e1 * t7, 0, t42, t13, t26 * t20 + t21 * t27, t51, t20 * t27 - t26 * t21, t1 * t9 + t5 * t2 - t22 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t38, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, -t38, 0, 0, 0, 0, 0, 0, -t11, 0, t12, -t51 - t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t10, 0, -t16, 0, 0, -t21 * t38, t20 * t38, 0, 0, t16, 0, t10, 0, 0, -t16, (-t20 * t8 - t21 * t5) * qJD(1), ((t6 - t35) * t21 + (t30 + t4) * t20) * qJD(1), qJD(4) * t49 + (-t20 * t5 + t21 * t8) * qJD(1), t3 * qJ(4) + t6 * qJD(4) - t5 * t8 + (-t48 + (-t4 - t39) * t20) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t19 * t24 - t23, t21 * t37 + (-t6 + t46) * qJD(3);];
tauc_reg = t18;
