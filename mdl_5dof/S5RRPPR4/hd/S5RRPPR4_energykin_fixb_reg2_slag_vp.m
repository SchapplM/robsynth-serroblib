% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:48
% EndTime: 2019-12-31 19:27:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (168->26), mult. (223->68), div. (0->0), fcn. (82->6), ass. (0->26)
t14 = qJD(1) + qJD(2);
t13 = t14 ^ 2;
t31 = t13 / 0.2e1;
t19 = sin(qJ(2));
t28 = pkin(1) * qJD(1);
t26 = t19 * t28;
t11 = t14 * qJ(3) + t26;
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t21 = cos(qJ(2));
t25 = t21 * t28;
t24 = qJD(3) - t25;
t8 = (-pkin(2) - pkin(3)) * t14 + t24;
t6 = t17 * t11 + t16 * t8;
t5 = -t16 * t11 + t17 * t8;
t3 = t14 * pkin(4) - t5;
t29 = t14 * t3;
t27 = qJD(5) * t14;
t22 = qJD(1) ^ 2;
t20 = cos(qJ(5));
t18 = sin(qJ(5));
t10 = -t14 * pkin(2) + t24;
t4 = -t14 * pkin(7) + t6;
t2 = t18 * qJD(4) + t20 * t4;
t1 = t20 * qJD(4) - t18 * t4;
t7 = [0, 0, 0, 0, 0, t22 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t14 * t25, -t14 * t26, 0, (t19 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t22, 0, 0, 0, t31, 0, 0, -t10 * t14, 0, t11 * t14, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t31, -t5 * t14, t6 * t14, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t18 ^ 2 * t31, t18 * t13 * t20, -t18 * t27, t20 ^ 2 * t31, -t20 * t27, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t20 * t29, -t2 * qJD(5) - t18 * t29, (t1 * t18 - t2 * t20) * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
