% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:04:39
% EndTime: 2021-01-15 17:04:40
% DurationCPUTime: 0.05s
% Computational Cost: add. (72->26), mult. (155->58), div. (0->0), fcn. (55->4), ass. (0->22)
t21 = qJD(1) ^ 2;
t31 = t21 / 0.2e1;
t18 = cos(pkin(7));
t25 = -pkin(1) * t18 - pkin(2);
t11 = qJD(3) + (-pkin(6) + t25) * qJD(1);
t19 = sin(qJ(4));
t20 = cos(qJ(4));
t30 = t20 * qJD(2) + t19 * t11;
t17 = sin(pkin(7));
t24 = -pkin(1) * t17 - qJ(3);
t10 = qJD(5) + (pkin(4) * t19 - t24) * qJD(1);
t29 = qJD(1) * t10;
t13 = t24 * qJD(1);
t28 = t13 * qJD(1);
t27 = qJ(5) * qJD(1);
t26 = qJD(1) * qJD(4);
t23 = -t19 * qJD(2) + t20 * t11;
t16 = qJD(2) ^ 2 / 0.2e1;
t12 = t25 * qJD(1) + qJD(3);
t7 = -t19 * t27 + t30;
t6 = qJD(4) * pkin(4) - t20 * t27 + t23;
t1 = [t31, 0, 0, t16 + (t17 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t21, t12 * qJD(1), -t28, t16 + t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t20 ^ 2 * t31, -t20 * t21 * t19, t20 * t26, -t19 * t26, qJD(4) ^ 2 / 0.2e1, t23 * qJD(4) - t19 * t28, -t30 * qJD(4) - t20 * t28, t6 * qJD(4) + t19 * t29, -t7 * qJD(4) + t20 * t29, (-t19 * t7 - t20 * t6) * qJD(1), t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t1;
