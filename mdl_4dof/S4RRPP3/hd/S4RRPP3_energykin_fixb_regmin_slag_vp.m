% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:36:03
% EndTime: 2021-01-15 10:36:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (101->25), mult. (269->60), div. (0->0), fcn. (143->4), ass. (0->24)
t29 = qJD(1) ^ 2;
t37 = t29 / 0.2e1;
t28 = cos(qJ(2));
t36 = t28 * t29;
t35 = pkin(5) + qJ(3);
t27 = sin(qJ(2));
t34 = qJD(1) * t27;
t21 = qJD(2) * pkin(2) - t35 * t34;
t33 = qJD(1) * t28;
t22 = t35 * t33;
t25 = sin(pkin(6));
t26 = cos(pkin(6));
t16 = t25 * t21 + t26 * t22;
t32 = qJD(1) * qJD(2);
t31 = t27 * t32;
t30 = t28 * t32;
t15 = t26 * t21 - t25 * t22;
t23 = qJD(3) + (-pkin(2) * t28 - pkin(1)) * qJD(1);
t19 = (t25 * t28 + t26 * t27) * qJD(1);
t18 = t25 * t34 - t26 * t33;
t14 = qJD(2) * qJ(4) + t16;
t13 = -qJD(2) * pkin(3) + qJD(4) - t15;
t12 = t18 * pkin(3) - t19 * qJ(4) + t23;
t1 = [t37, 0, 0, t27 ^ 2 * t37, t27 * t36, t31, t30, qJD(2) ^ 2 / 0.2e1, pkin(1) * t36 - pkin(5) * t31, -t29 * pkin(1) * t27 - pkin(5) * t30, t15 * qJD(2) + t23 * t18, -t16 * qJD(2) + t23 * t19, -t15 * t19 - t16 * t18, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, -t13 * qJD(2) + t12 * t18, t13 * t19 - t14 * t18, t14 * qJD(2) - t12 * t19, t14 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg = t1;
