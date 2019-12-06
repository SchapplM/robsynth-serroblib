% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:40
% EndTime: 2019-12-05 15:33:40
% DurationCPUTime: 0.14s
% Computational Cost: add. (104->29), mult. (265->74), div. (0->0), fcn. (145->6), ass. (0->33)
t30 = qJD(2) ^ 2;
t23 = t30 / 0.2e1;
t26 = sin(qJ(4));
t28 = cos(qJ(4));
t29 = cos(qJ(2));
t12 = qJD(2) * pkin(2) + t29 * qJD(1);
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t27 = sin(qJ(2));
t37 = qJD(1) * t27;
t10 = t24 * t12 + t25 * t37;
t8 = qJD(2) * pkin(6) + t10;
t4 = t26 * qJD(3) + t28 * t8;
t36 = qJD(2) * t26;
t35 = qJD(2) * t28;
t34 = qJ(5) * qJD(2);
t33 = qJD(1) * qJD(2);
t32 = qJD(2) * qJD(4);
t9 = t25 * t12 - t24 * t37;
t31 = qJD(1) ^ 2;
t22 = qJD(4) ^ 2 / 0.2e1;
t21 = t28 * qJD(3);
t19 = t28 * t32;
t18 = t26 * t32;
t17 = t28 ^ 2 * t23;
t16 = t26 ^ 2 * t23;
t15 = t26 * t30 * t28;
t7 = -qJD(2) * pkin(3) - t9;
t5 = qJD(5) + (-pkin(4) * t28 - pkin(3)) * qJD(2) - t9;
t3 = -t26 * t8 + t21;
t2 = t28 * t34 + t4;
t1 = qJD(4) * pkin(4) + t21 + (-t8 - t34) * t26;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t31 / 0.2e1, 0, 0, 0, 0, 0, t23, t29 * t33, -t27 * t33, 0, (t27 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1) * t31, 0, 0, 0, 0, 0, t23, t9 * qJD(2), -t10 * qJD(2), 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t16, t15, t18, t17, t19, t22, t3 * qJD(4) - t35 * t7, -t4 * qJD(4) + t36 * t7, (-t26 * t3 + t28 * t4) * qJD(2), t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t16, t15, t18, t17, t19, t22, t1 * qJD(4) - t35 * t5, -t2 * qJD(4) + t36 * t5, (-t1 * t26 + t2 * t28) * qJD(2), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
