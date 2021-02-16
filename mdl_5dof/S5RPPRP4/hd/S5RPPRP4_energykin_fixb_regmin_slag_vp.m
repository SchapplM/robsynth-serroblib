% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:13
% EndTime: 2021-01-15 17:13:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (96->28), mult. (179->60), div. (0->0), fcn. (66->4), ass. (0->23)
t27 = qJD(1) ^ 2;
t34 = t27 / 0.2e1;
t16 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t22 = sin(pkin(7));
t23 = cos(pkin(7));
t30 = qJ(2) * qJD(1);
t14 = t22 * t16 + t23 * t30;
t12 = -qJD(1) * pkin(6) + t14;
t25 = sin(qJ(4));
t26 = cos(qJ(4));
t33 = t25 * qJD(3) + t26 * t12;
t32 = qJD(1) * t25;
t31 = qJD(1) * t26;
t29 = qJ(5) * qJD(1);
t28 = qJD(1) * qJD(4);
t13 = t23 * t16 - t22 * t30;
t11 = qJD(1) * pkin(3) - t13;
t21 = t26 * qJD(3);
t19 = -qJD(1) * pkin(1) + qJD(2);
t9 = pkin(4) * t31 + qJD(5) + t11;
t8 = -t26 * t29 + t33;
t7 = qJD(4) * pkin(4) + t21 + (-t12 + t29) * t25;
t1 = [t34, 0, 0, -t19 * qJD(1), t27 * qJ(2), qJ(2) ^ 2 * t34 + t19 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t25 ^ 2 * t34, t25 * t27 * t26, -t25 * t28, -t26 * t28, qJD(4) ^ 2 / 0.2e1, t11 * t31 + (-t25 * t12 + t21) * qJD(4), -t33 * qJD(4) - t11 * t32, t7 * qJD(4) + t9 * t31, -t8 * qJD(4) - t9 * t32, (t25 * t7 - t26 * t8) * qJD(1), t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg = t1;
