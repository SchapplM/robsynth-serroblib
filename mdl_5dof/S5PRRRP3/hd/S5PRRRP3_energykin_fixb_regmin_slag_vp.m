% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:25
% EndTime: 2021-01-15 16:23:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (109->28), mult. (265->65), div. (0->0), fcn. (160->4), ass. (0->24)
t29 = qJD(2) ^ 2;
t38 = t29 / 0.2e1;
t37 = cos(qJ(4));
t28 = cos(qJ(3));
t36 = t28 * t29;
t24 = t28 * qJD(1);
t27 = sin(qJ(3));
t33 = qJD(2) * t27;
t16 = qJD(3) * pkin(3) + t24 + (-pkin(7) - pkin(6)) * t33;
t32 = qJD(2) * t28;
t34 = pkin(6) * t32 + t27 * qJD(1);
t17 = pkin(7) * t32 + t34;
t26 = sin(qJ(4));
t35 = t26 * t16 + t37 * t17;
t31 = qJD(2) * qJD(3);
t30 = t37 * t16 - t26 * t17;
t20 = (-pkin(3) * t28 - pkin(2)) * qJD(2);
t25 = qJD(3) + qJD(4);
t19 = (t26 * t28 + t37 * t27) * qJD(2);
t18 = t26 * t33 - t37 * t32;
t12 = t18 * pkin(4) + qJD(5) + t20;
t11 = -t18 * qJ(5) + t35;
t10 = t25 * pkin(4) - t19 * qJ(5) + t30;
t1 = [qJD(1) ^ 2 / 0.2e1, t38, 0, 0, t27 ^ 2 * t38, t27 * t36, t27 * t31, t28 * t31, qJD(3) ^ 2 / 0.2e1, pkin(2) * t36 + (-pkin(6) * t33 + t24) * qJD(3), -t29 * pkin(2) * t27 - t34 * qJD(3), t19 ^ 2 / 0.2e1, -t19 * t18, t19 * t25, -t18 * t25, t25 ^ 2 / 0.2e1, t20 * t18 + t30 * t25, t20 * t19 - t35 * t25, t10 * t25 + t12 * t18, -t11 * t25 + t12 * t19, -t10 * t19 - t11 * t18, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1;];
T_reg = t1;
