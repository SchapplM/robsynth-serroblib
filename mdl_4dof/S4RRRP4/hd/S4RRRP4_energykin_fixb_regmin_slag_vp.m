% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:21
% EndTime: 2021-01-15 14:30:21
% DurationCPUTime: 0.06s
% Computational Cost: add. (91->24), mult. (245->59), div. (0->0), fcn. (142->4), ass. (0->25)
t25 = qJD(1) ^ 2;
t36 = t25 / 0.2e1;
t35 = pkin(6) + pkin(5);
t34 = cos(qJ(3));
t24 = cos(qJ(2));
t33 = t24 * t25;
t23 = sin(qJ(2));
t31 = qJD(1) * t23;
t17 = qJD(2) * pkin(2) - t35 * t31;
t30 = qJD(1) * t24;
t18 = t35 * t30;
t22 = sin(qJ(3));
t32 = t22 * t17 + t34 * t18;
t29 = qJD(1) * qJD(2);
t28 = t23 * t29;
t27 = t24 * t29;
t26 = t34 * t17 - t22 * t18;
t19 = (-pkin(2) * t24 - pkin(1)) * qJD(1);
t21 = qJD(2) + qJD(3);
t15 = (t22 * t24 + t34 * t23) * qJD(1);
t14 = t22 * t31 - t34 * t30;
t11 = t14 * pkin(3) + qJD(4) + t19;
t10 = -t14 * qJ(4) + t32;
t9 = t21 * pkin(3) - t15 * qJ(4) + t26;
t1 = [t36, 0, 0, t23 ^ 2 * t36, t23 * t33, t28, t27, qJD(2) ^ 2 / 0.2e1, pkin(1) * t33 - pkin(5) * t28, -t25 * pkin(1) * t23 - pkin(5) * t27, t15 ^ 2 / 0.2e1, -t15 * t14, t15 * t21, -t14 * t21, t21 ^ 2 / 0.2e1, t19 * t14 + t26 * t21, t19 * t15 - t32 * t21, t11 * t14 + t9 * t21, -t10 * t21 + t11 * t15, -t10 * t14 - t9 * t15, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1;];
T_reg = t1;
