% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPR8
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:41
% EndTime: 2019-12-31 17:42:41
% DurationCPUTime: 0.10s
% Computational Cost: add. (150->23), mult. (277->71), div. (0->0), fcn. (172->8), ass. (0->27)
t16 = qJD(2) + qJD(3);
t15 = t16 ^ 2;
t14 = t15 / 0.2e1;
t24 = cos(qJ(2));
t13 = qJD(2) * pkin(2) + t24 * qJD(1);
t20 = sin(qJ(3));
t23 = cos(qJ(3));
t21 = sin(qJ(2));
t28 = qJD(1) * t21;
t11 = t20 * t13 + t23 * t28;
t17 = sin(pkin(9));
t18 = cos(pkin(9));
t10 = t23 * t13 - t20 * t28;
t8 = t16 * pkin(3) + t10;
t6 = t18 * t11 + t17 * t8;
t5 = -t17 * t11 + t18 * t8;
t3 = -t16 * pkin(4) - t5;
t29 = t16 * t3;
t27 = qJD(5) * t16;
t26 = qJD(1) * qJD(2);
t25 = qJD(1) ^ 2;
t22 = cos(qJ(5));
t19 = sin(qJ(5));
t4 = t16 * pkin(7) + t6;
t2 = t19 * qJD(4) + t22 * t4;
t1 = t22 * qJD(4) - t19 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t25 / 0.2e1, 0, 0, 0, 0, 0, qJD(2) ^ 2 / 0.2e1, t24 * t26, -t21 * t26, 0, (t21 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1) * t25, 0, 0, 0, 0, 0, t14, t10 * t16, -t11 * t16, 0, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t14, t5 * t16, -t6 * t16, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t19 ^ 2 * t14, t19 * t15 * t22, t19 * t27, t22 ^ 2 * t14, t22 * t27, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t22 * t29, -t2 * qJD(5) + t19 * t29, (-t1 * t19 + t2 * t22) * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
