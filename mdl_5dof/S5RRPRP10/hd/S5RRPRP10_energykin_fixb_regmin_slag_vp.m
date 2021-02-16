% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:44
% EndTime: 2021-01-15 21:02:44
% DurationCPUTime: 0.08s
% Computational Cost: add. (148->36), mult. (326->78), div. (0->0), fcn. (161->4), ass. (0->30)
t31 = qJD(1) ^ 2;
t43 = t31 / 0.2e1;
t42 = -pkin(7) - pkin(2);
t30 = cos(qJ(2));
t41 = t30 * t31;
t28 = sin(qJ(2));
t33 = -qJ(3) * t28 - pkin(1);
t14 = (t42 * t30 + t33) * qJD(1);
t38 = t28 * qJD(1);
t37 = pkin(6) * t38 + qJD(3);
t15 = pkin(3) * t38 + t42 * qJD(2) + t37;
t27 = sin(qJ(4));
t29 = cos(qJ(4));
t40 = t29 * t14 + t27 * t15;
t39 = qJD(1) * t30;
t21 = -pkin(6) * t39 - qJD(2) * qJ(3);
t36 = qJD(1) * qJD(2);
t16 = pkin(3) * t39 - t21;
t35 = t28 * t36;
t34 = t30 * t36;
t32 = -t27 * t14 + t29 * t15;
t22 = qJD(4) + t38;
t20 = -qJD(2) * pkin(2) + t37;
t19 = t29 * qJD(2) - t27 * t39;
t18 = t27 * qJD(2) + t29 * t39;
t17 = (-pkin(2) * t30 + t33) * qJD(1);
t10 = t18 * pkin(4) + qJD(5) + t16;
t9 = -t18 * qJ(5) + t40;
t8 = t22 * pkin(4) - t19 * qJ(5) + t32;
t1 = [t43, 0, 0, t28 ^ 2 * t43, t28 * t41, t35, t34, qJD(2) ^ 2 / 0.2e1, pkin(1) * t41 - pkin(6) * t35, -t31 * pkin(1) * t28 - pkin(6) * t34, (t20 * t28 - t21 * t30) * qJD(1), t20 * qJD(2) + t17 * t39, -t21 * qJD(2) - t17 * t38, t17 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t18, t19 * t22, -t18 * t22, t22 ^ 2 / 0.2e1, t16 * t18 + t32 * t22, t16 * t19 - t40 * t22, t10 * t18 + t8 * t22, t10 * t19 - t9 * t22, -t9 * t18 - t8 * t19, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t1;
