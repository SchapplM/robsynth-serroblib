% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:49
% EndTime: 2019-12-05 16:15:49
% DurationCPUTime: 0.11s
% Computational Cost: add. (154->29), mult. (254->77), div. (0->0), fcn. (133->6), ass. (0->30)
t18 = qJD(2) + qJD(3);
t17 = t18 ^ 2;
t34 = t17 / 0.2e1;
t33 = cos(qJ(5));
t20 = sin(pkin(9));
t32 = t18 * t20;
t21 = cos(pkin(9));
t31 = t18 * t21;
t23 = sin(qJ(3));
t30 = pkin(2) * qJD(2);
t29 = t23 * t30;
t13 = t18 * qJ(4) + t29;
t6 = t20 * qJD(1) + t21 * t13;
t24 = cos(qJ(3));
t28 = t24 * t30;
t27 = qJD(4) - t28;
t25 = qJD(2) ^ 2;
t22 = sin(qJ(5));
t19 = qJD(1) ^ 2 / 0.2e1;
t16 = t21 * qJD(1);
t12 = -t18 * pkin(3) + t27;
t10 = (t33 * t20 + t21 * t22) * t18;
t8 = t22 * t32 - t33 * t31;
t7 = (-pkin(4) * t21 - pkin(3)) * t18 + t27;
t5 = -t20 * t13 + t16;
t4 = pkin(7) * t31 + t6;
t3 = t16 + (-pkin(7) * t18 - t13) * t20;
t2 = t22 * t3 + t33 * t4;
t1 = -t22 * t4 + t33 * t3;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, t25 / 0.2e1, 0, 0, 0, t19, 0, 0, 0, 0, 0, t34, t18 * t28, -t18 * t29, 0, t19 + (t23 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t25, t20 ^ 2 * t34, t20 * t17 * t21, 0, t21 ^ 2 * t34, 0, 0, -t12 * t31, t12 * t32, (-t20 * t5 + t21 * t6) * t18, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * qJD(5), t8 ^ 2 / 0.2e1, -t8 * qJD(5), qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t7 * t8, -t2 * qJD(5) + t7 * t10, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
