% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:25
% EndTime: 2019-03-08 18:26:25
% DurationCPUTime: 0.13s
% Computational Cost: add. (91->38), mult. (374->81), div. (0->0), fcn. (210->4), ass. (0->35)
t23 = qJD(1) ^ 2;
t37 = t23 / 0.2e1;
t20 = sin(pkin(4));
t36 = t20 ^ 2 * t23;
t19 = sin(pkin(6));
t35 = t19 * t20;
t21 = cos(pkin(6));
t34 = t20 * t21;
t33 = qJD(1) * t20;
t26 = qJ(2) * t33;
t22 = cos(pkin(4));
t32 = qJD(1) * t22;
t29 = pkin(1) * t32;
t8 = t19 * t29 + t21 * t26;
t13 = t19 * t26;
t31 = qJD(3) + t13;
t30 = t20 * t22 * t23;
t28 = t36 / 0.2e1;
t27 = -pkin(1) * t21 - pkin(2);
t25 = -qJ(3) * t19 - pkin(1);
t24 = t21 * t30;
t17 = t22 ^ 2 * t37;
t16 = -pkin(1) * t33 + qJD(2);
t12 = t21 ^ 2 * t28;
t11 = t19 ^ 2 * t28;
t10 = t19 * t30;
t9 = t19 * t21 * t36;
t7 = t21 * t29 - t13;
t6 = -qJ(3) * t32 - t8;
t5 = qJD(2) + (-pkin(2) * t21 + t25) * t33;
t4 = t27 * t32 + t31;
t3 = qJD(2) + ((-pkin(2) - qJ(4)) * t21 + t25) * t33;
t2 = qJD(4) + (pkin(3) * t34 + qJ(3) * t22) * qJD(1) + t8;
t1 = (pkin(3) * t35 + (-qJ(4) + t27) * t22) * qJD(1) + t31;
t14 = [0, 0, 0, 0, 0, t37, 0, 0, 0, 0, t11, t9, t10, t12, t24, t17 (-t16 * t34 + t22 * t7) * qJD(1) (t16 * t35 - t22 * t8) * qJD(1) (-t19 * t7 + t21 * t8) * t33, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t17, -t10, -t24, t11, t9, t12 (t19 * t4 - t21 * t6) * t33 (t22 * t4 + t5 * t34) * qJD(1) (-t22 * t6 - t5 * t35) * qJD(1), t5 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t17, -t24, t10, t12, -t9, t11 (t1 * t19 + t2 * t21) * t33 (t2 * t22 - t3 * t35) * qJD(1) (-t1 * t22 - t3 * t34) * qJD(1), t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg  = t14;
