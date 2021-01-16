% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:35
% EndTime: 2021-01-15 10:20:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (45->19), mult. (127->53), div. (0->0), fcn. (51->4), ass. (0->19)
t19 = qJD(1) ^ 2;
t27 = t19 / 0.2e1;
t15 = sin(pkin(6));
t11 = (pkin(1) * t15 + pkin(5)) * qJD(1);
t17 = sin(qJ(3));
t18 = cos(qJ(3));
t26 = t17 * qJD(2) + t18 * t11;
t25 = qJD(1) * t17;
t24 = qJD(1) * t18;
t23 = qJ(4) * qJD(1);
t22 = qJD(1) * qJD(3);
t16 = cos(pkin(6));
t21 = -pkin(1) * t16 - pkin(2);
t14 = t18 * qJD(2);
t12 = t21 * qJD(1);
t9 = qJD(4) + (-pkin(3) * t18 + t21) * qJD(1);
t8 = t18 * t23 + t26;
t7 = qJD(3) * pkin(3) + t14 + (-t11 - t23) * t17;
t1 = [t27, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t15 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t19, t17 ^ 2 * t27, t17 * t19 * t18, t17 * t22, t18 * t22, qJD(3) ^ 2 / 0.2e1, -t12 * t24 + (-t17 * t11 + t14) * qJD(3), -t26 * qJD(3) + t12 * t25, t7 * qJD(3) - t9 * t24, -t8 * qJD(3) + t9 * t25, (-t17 * t7 + t18 * t8) * qJD(1), t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg = t1;
