% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:04
% EndTime: 2019-12-31 20:11:04
% DurationCPUTime: 0.14s
% Computational Cost: add. (202->46), mult. (485->95), div. (0->0), fcn. (232->4), ass. (0->42)
t40 = qJD(1) ^ 2;
t51 = t40 / 0.2e1;
t50 = -pkin(2) - pkin(7);
t39 = cos(qJ(2));
t37 = sin(qJ(2));
t42 = -qJ(3) * t37 - pkin(1);
t10 = (t50 * t39 + t42) * qJD(1);
t47 = t37 * qJD(1);
t46 = pkin(6) * t47 + qJD(3);
t13 = pkin(3) * t47 + t50 * qJD(2) + t46;
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t4 = t38 * t10 + t36 * t13;
t49 = t39 * t40;
t48 = qJD(1) * t39;
t22 = -pkin(6) * t48 - qJD(2) * qJ(3);
t45 = qJD(1) * qJD(2);
t14 = pkin(3) * t48 - t22;
t3 = -t36 * t10 + t38 * t13;
t44 = t37 * t45;
t43 = t39 * t45;
t35 = t39 ^ 2;
t34 = t37 ^ 2;
t32 = qJD(2) ^ 2 / 0.2e1;
t28 = t35 * t51;
t27 = t34 * t51;
t26 = qJD(4) + t47;
t25 = t37 * t49;
t23 = t26 ^ 2 / 0.2e1;
t21 = -qJD(2) * pkin(2) + t46;
t20 = t38 * qJD(2) - t36 * t48;
t18 = t36 * qJD(2) + t38 * t48;
t17 = t20 ^ 2 / 0.2e1;
t16 = t18 ^ 2 / 0.2e1;
t15 = (-pkin(2) * t39 + t42) * qJD(1);
t12 = t20 * t26;
t11 = t18 * t26;
t6 = t20 * t18;
t5 = t18 * pkin(4) + qJD(5) + t14;
t2 = -t18 * qJ(5) + t4;
t1 = t26 * pkin(4) - t20 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, t51, 0, 0, 0, 0, t27, t25, t44, t28, t43, t32, pkin(1) * t49 - pkin(6) * t44, -t40 * pkin(1) * t37 - pkin(6) * t43, (t34 + t35) * t40 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t35 / 0.2e1 + t34 / 0.2e1) * pkin(6) ^ 2) * t40, t32, -t44, -t43, t27, t25, t28, (t21 * t37 - t22 * t39) * qJD(1), t21 * qJD(2) + t15 * t48, -t22 * qJD(2) - t15 * t47, t15 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t17, -t6, t12, t16, -t11, t23, t14 * t18 + t3 * t26, t14 * t20 - t4 * t26, -t4 * t18 - t3 * t20, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17, -t6, t12, t16, -t11, t23, t1 * t26 + t5 * t18, -t2 * t26 + t5 * t20, -t1 * t20 - t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
