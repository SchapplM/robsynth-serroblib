% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:57
% EndTime: 2019-12-05 17:07:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (168->31), mult. (278->85), div. (0->0), fcn. (141->6), ass. (0->31)
t19 = qJD(2) + qJD(3);
t17 = t19 ^ 2;
t35 = t17 / 0.2e1;
t34 = cos(qJ(5));
t22 = sin(qJ(4));
t33 = t19 * t22;
t24 = cos(qJ(4));
t32 = t19 * t24;
t23 = sin(qJ(3));
t31 = pkin(2) * qJD(2);
t29 = t23 * t31;
t12 = t19 * pkin(7) + t29;
t6 = t22 * qJD(1) + t24 * t12;
t30 = qJD(4) * t19;
t25 = cos(qJ(3));
t28 = t25 * t31;
t26 = qJD(2) ^ 2;
t21 = sin(qJ(5));
t20 = qJD(1) ^ 2 / 0.2e1;
t18 = qJD(4) + qJD(5);
t16 = t24 * qJD(1);
t13 = -t19 * pkin(3) - t28;
t10 = -t28 + (-pkin(4) * t24 - pkin(3)) * t19;
t9 = (t21 * t24 + t34 * t22) * t19;
t7 = t21 * t33 - t34 * t32;
t5 = -t22 * t12 + t16;
t4 = pkin(8) * t32 + t6;
t3 = qJD(4) * pkin(4) + t16 + (-pkin(8) * t19 - t12) * t22;
t2 = t21 * t3 + t34 * t4;
t1 = -t21 * t4 + t34 * t3;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, t26 / 0.2e1, 0, 0, 0, t20, 0, 0, 0, 0, 0, t35, t19 * t28, -t19 * t29, 0, t20 + (t23 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t26, t22 ^ 2 * t35, t22 * t17 * t24, t22 * t30, t24 ^ 2 * t35, t24 * t30, qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) - t13 * t32, -t6 * qJD(4) + t13 * t33, (-t22 * t5 + t24 * t6) * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t18, t7 ^ 2 / 0.2e1, -t7 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t10 * t7, t10 * t9 - t2 * t18, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
