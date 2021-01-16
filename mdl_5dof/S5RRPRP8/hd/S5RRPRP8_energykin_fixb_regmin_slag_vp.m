% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP8
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
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:18
% EndTime: 2021-01-15 20:52:18
% DurationCPUTime: 0.08s
% Computational Cost: add. (145->36), mult. (330->79), div. (0->0), fcn. (168->4), ass. (0->28)
t33 = qJD(1) ^ 2;
t43 = t33 / 0.2e1;
t30 = sin(qJ(2));
t40 = qJD(1) * t30;
t38 = pkin(6) * t40 + qJD(3);
t12 = -pkin(7) * t40 + (-pkin(2) - pkin(3)) * qJD(2) + t38;
t32 = cos(qJ(2));
t39 = qJD(1) * t32;
t19 = pkin(6) * t39 + qJD(2) * qJ(3);
t16 = -pkin(7) * t39 + t19;
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t42 = t29 * t12 + t31 * t16;
t41 = t32 * t33;
t37 = qJD(1) * qJD(2);
t17 = -qJD(1) * pkin(1) - pkin(2) * t39 - qJ(3) * t40;
t36 = t30 * t37;
t35 = t32 * t37;
t34 = t31 * t12 - t29 * t16;
t11 = pkin(3) * t39 - t17;
t25 = qJD(2) - qJD(4);
t18 = -qJD(2) * pkin(2) + t38;
t15 = (-t29 * t32 + t30 * t31) * qJD(1);
t14 = (t29 * t30 + t31 * t32) * qJD(1);
t8 = t14 * pkin(4) + qJD(5) + t11;
t7 = -t14 * qJ(5) + t42;
t6 = -t25 * pkin(4) - t15 * qJ(5) + t34;
t1 = [t43, 0, 0, t30 ^ 2 * t43, t30 * t41, t36, t35, qJD(2) ^ 2 / 0.2e1, pkin(1) * t41 - pkin(6) * t36, -t33 * pkin(1) * t30 - pkin(6) * t35, -t18 * qJD(2) - t17 * t39, (t18 * t30 + t19 * t32) * qJD(1), t19 * qJD(2) - t17 * t40, t19 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t14, -t15 * t25, t14 * t25, t25 ^ 2 / 0.2e1, t11 * t14 - t34 * t25, t11 * t15 + t42 * t25, t8 * t14 - t6 * t25, t8 * t15 + t7 * t25, -t7 * t14 - t6 * t15, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg = t1;
