% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:30
% EndTime: 2019-12-31 20:01:30
% DurationCPUTime: 0.17s
% Computational Cost: add. (333->46), mult. (861->107), div. (0->0), fcn. (553->6), ass. (0->43)
t38 = qJD(1) ^ 2;
t52 = t38 / 0.2e1;
t36 = sin(qJ(2));
t45 = qJD(1) * t36;
t47 = pkin(6) + qJ(3);
t26 = qJD(2) * pkin(2) - t47 * t45;
t37 = cos(qJ(2));
t44 = qJD(1) * t37;
t27 = t47 * t44;
t34 = sin(pkin(8));
t46 = cos(pkin(8));
t13 = t34 * t26 + t46 * t27;
t11 = qJD(2) * pkin(7) + t13;
t35 = sin(qJ(4));
t51 = cos(qJ(4));
t21 = t34 * t45 - t46 * t44;
t23 = (t34 * t37 + t46 * t36) * qJD(1);
t28 = qJD(3) + (-pkin(2) * t37 - pkin(1)) * qJD(1);
t7 = t21 * pkin(3) - t23 * pkin(7) + t28;
t4 = t51 * t11 + t35 * t7;
t15 = -t51 * qJD(2) + t35 * t23;
t17 = t35 * qJD(2) + t51 * t23;
t50 = t17 * t15;
t20 = qJD(4) + t21;
t49 = t20 * t15;
t48 = t37 * t38;
t43 = t15 ^ 2 / 0.2e1;
t42 = qJD(1) * qJD(2);
t41 = t36 * t42;
t40 = t37 * t42;
t12 = t46 * t26 - t34 * t27;
t3 = -t35 * t11 + t51 * t7;
t10 = -qJD(2) * pkin(3) - t12;
t33 = t37 ^ 2;
t32 = t36 ^ 2;
t31 = qJD(2) ^ 2 / 0.2e1;
t18 = t20 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t8 = t17 * t20;
t5 = t15 * pkin(4) - t17 * qJ(5) + t10;
t2 = t20 * qJ(5) + t4;
t1 = -t20 * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t52, 0, 0, 0, 0, t32 * t52, t36 * t48, t41, t33 * t52, t40, t31, pkin(1) * t48 - pkin(6) * t41, -t38 * pkin(1) * t36 - pkin(6) * t40, (t32 + t33) * t38 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t33 / 0.2e1 + t32 / 0.2e1) * pkin(6) ^ 2) * t38, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * qJD(2), t21 ^ 2 / 0.2e1, -t21 * qJD(2), t31, t12 * qJD(2) + t28 * t21, -t13 * qJD(2) + t28 * t23, -t12 * t23 - t13 * t21, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t14, -t50, t8, t43, -t49, t18, t10 * t15 + t3 * t20, t10 * t17 - t4 * t20, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t14, t8, t50, t18, t49, t43, -t1 * t20 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t20, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
