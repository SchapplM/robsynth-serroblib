% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRP6
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
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:13
% EndTime: 2021-01-15 14:39:13
% DurationCPUTime: 0.05s
% Computational Cost: add. (91->25), mult. (227->60), div. (0->0), fcn. (124->4), ass. (0->24)
t29 = qJD(1) ^ 2;
t39 = t29 / 0.2e1;
t38 = cos(qJ(3));
t28 = cos(qJ(2));
t37 = t28 * t29;
t27 = sin(qJ(2));
t16 = (-pkin(2) * t28 - pkin(6) * t27 - pkin(1)) * qJD(1);
t34 = t28 * qJD(1);
t21 = pkin(5) * t34 + qJD(2) * pkin(6);
t26 = sin(qJ(3));
t36 = t26 * t16 + t38 * t21;
t35 = qJD(1) * t27;
t33 = qJD(1) * qJD(2);
t32 = t27 * t33;
t31 = t28 * t33;
t30 = t38 * t16 - t26 * t21;
t20 = -qJD(2) * pkin(2) + pkin(5) * t35;
t22 = -qJD(3) + t34;
t18 = t26 * qJD(2) + t38 * t35;
t17 = -t38 * qJD(2) + t26 * t35;
t13 = t17 * pkin(3) + qJD(4) + t20;
t12 = -t17 * qJ(4) + t36;
t11 = -t22 * pkin(3) - t18 * qJ(4) + t30;
t1 = [t39, 0, 0, t27 ^ 2 * t39, t27 * t37, t32, t31, qJD(2) ^ 2 / 0.2e1, pkin(1) * t37 - pkin(5) * t32, -t29 * pkin(1) * t27 - pkin(5) * t31, t18 ^ 2 / 0.2e1, -t18 * t17, -t18 * t22, t17 * t22, t22 ^ 2 / 0.2e1, t20 * t17 - t30 * t22, t20 * t18 + t36 * t22, -t11 * t22 + t13 * t17, t12 * t22 + t13 * t18, -t11 * t18 - t12 * t17, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg = t1;
