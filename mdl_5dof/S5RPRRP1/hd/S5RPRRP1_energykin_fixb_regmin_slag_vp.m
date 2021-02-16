% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:18
% EndTime: 2021-01-15 12:27:18
% DurationCPUTime: 0.07s
% Computational Cost: add. (132->30), mult. (274->67), div. (0->0), fcn. (142->4), ass. (0->24)
t27 = qJD(1) ^ 2;
t34 = t27 / 0.2e1;
t26 = cos(qJ(3));
t18 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t28 = -pkin(7) * qJD(1) + t18;
t13 = qJD(3) * pkin(3) + t28 * t26;
t24 = sin(qJ(3));
t14 = t28 * t24;
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t33 = t23 * t13 + t25 * t14;
t17 = (pkin(3) * t24 + qJ(2)) * qJD(1);
t32 = t27 * qJ(2);
t31 = qJD(3) * t18;
t30 = qJD(1) * qJD(3);
t29 = t25 * t13 - t23 * t14;
t21 = qJD(3) + qJD(4);
t20 = -qJD(1) * pkin(1) + qJD(2);
t16 = (-t23 * t24 + t25 * t26) * qJD(1);
t15 = (t23 * t26 + t24 * t25) * qJD(1);
t9 = t15 * pkin(4) + qJD(5) + t17;
t8 = -t15 * qJ(5) + t33;
t7 = t21 * pkin(4) - t16 * qJ(5) + t29;
t1 = [t34, 0, 0, t20 * qJD(1), t32, qJ(2) ^ 2 * t34 + t20 ^ 2 / 0.2e1, t26 ^ 2 * t34, -t26 * t27 * t24, t26 * t30, -t24 * t30, qJD(3) ^ 2 / 0.2e1, t24 * t32 + t26 * t31, -t24 * t31 + t26 * t32, t16 ^ 2 / 0.2e1, -t15 * t16, t21 * t16, -t15 * t21, t21 ^ 2 / 0.2e1, t17 * t15 + t29 * t21, t17 * t16 - t33 * t21, t9 * t15 + t7 * t21, t9 * t16 - t8 * t21, -t8 * t15 - t7 * t16, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg = t1;
