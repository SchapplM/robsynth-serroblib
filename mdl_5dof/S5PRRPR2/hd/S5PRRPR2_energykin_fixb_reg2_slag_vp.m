% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPR2
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:34
% EndTime: 2019-12-05 16:17:34
% DurationCPUTime: 0.11s
% Computational Cost: add. (139->26), mult. (231->71), div. (0->0), fcn. (108->6), ass. (0->31)
t14 = qJD(2) + qJD(3);
t12 = t14 ^ 2;
t35 = t12 / 0.2e1;
t16 = sin(pkin(9));
t34 = t16 ^ 2 * t12;
t33 = t14 * t16;
t17 = cos(pkin(9));
t32 = t17 * t14;
t31 = pkin(2) * qJD(2);
t19 = sin(qJ(3));
t27 = t19 * t31;
t8 = t14 * qJ(4) + t27;
t4 = -t17 * qJD(1) + t16 * t8;
t30 = t4 ^ 2 / 0.2e1;
t18 = sin(qJ(5));
t29 = t18 * t33;
t20 = cos(qJ(5));
t28 = t20 * t33;
t21 = cos(qJ(3));
t26 = t21 * t31;
t25 = t34 / 0.2e1;
t24 = qJD(4) - t26;
t22 = qJD(2) ^ 2;
t15 = qJD(1) ^ 2 / 0.2e1;
t9 = -qJD(5) + t32;
t7 = -t14 * pkin(3) + t24;
t6 = t16 * qJD(1) + t17 * t8;
t3 = (-pkin(4) * t17 - pkin(7) * t16 - pkin(3)) * t14 + t24;
t2 = t18 * t3 + t20 * t6;
t1 = -t18 * t6 + t20 * t3;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, t22 / 0.2e1, 0, 0, 0, t15, 0, 0, 0, 0, 0, t35, t14 * t26, -t14 * t27, 0, t15 + (t19 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t22, t25, t16 * t12 * t17, 0, t17 ^ 2 * t35, 0, 0, -t7 * t32, t7 * t33, (t16 * t4 + t17 * t6) * t14, t6 ^ 2 / 0.2e1 + t30 + t7 ^ 2 / 0.2e1, t20 ^ 2 * t25, -t20 * t18 * t34, -t9 * t28, t18 ^ 2 * t25, t9 * t29, t9 ^ 2 / 0.2e1, -t1 * t9 + t4 * t29, t2 * t9 + t4 * t28, (-t1 * t20 - t18 * t2) * t33, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t30;];
T_reg = t5;
