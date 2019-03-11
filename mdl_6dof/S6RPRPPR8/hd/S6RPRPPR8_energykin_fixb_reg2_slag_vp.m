% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:59:58
% EndTime: 2019-03-09 02:59:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (222->54), mult. (440->104), div. (0->0), fcn. (162->4), ass. (0->43)
t34 = qJD(1) ^ 2;
t26 = t34 / 0.2e1;
t46 = -pkin(3) - pkin(4);
t17 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t31 = sin(qJ(3));
t11 = qJD(3) * qJ(4) + t31 * t17;
t45 = t34 * qJ(2);
t44 = qJD(1) * t31;
t43 = qJD(3) * t17;
t33 = cos(qJ(3));
t42 = t33 * qJD(1);
t22 = qJ(4) * t42;
t41 = qJD(5) + t22;
t40 = qJ(5) * qJD(1);
t39 = -pkin(8) + t46;
t38 = qJD(1) * qJD(3);
t37 = t33 * t34 * t31;
t36 = t31 * t38;
t8 = -t31 * t40 - t11;
t35 = qJD(4) + (-t17 - t40) * t33;
t32 = cos(qJ(6));
t30 = sin(qJ(6));
t29 = t33 ^ 2;
t28 = t31 ^ 2;
t25 = qJD(3) ^ 2 / 0.2e1;
t24 = qJ(2) ^ 2 * t26;
t23 = -pkin(1) * qJD(1) + qJD(2);
t21 = t33 * t38;
t20 = t29 * t26;
t19 = t28 * t26;
t18 = qJD(6) + t42;
t14 = -t30 * qJD(3) + t32 * t44;
t12 = t32 * qJD(3) + t30 * t44;
t10 = -t22 + (pkin(3) * t31 + qJ(2)) * qJD(1);
t9 = -qJD(3) * pkin(3) - t33 * t17 + qJD(4);
t7 = (t46 * t31 - qJ(2)) * qJD(1) + t41;
t6 = qJD(3) * pkin(5) - t8;
t5 = t46 * qJD(3) + t35;
t4 = t39 * qJD(3) + t35;
t3 = (pkin(5) * t33 + t39 * t31 - qJ(2)) * qJD(1) + t41;
t2 = t30 * t3 + t32 * t4;
t1 = t32 * t3 - t30 * t4;
t13 = [0, 0, 0, 0, 0, t26, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, t23 * qJD(1), t45, t24 + t23 ^ 2 / 0.2e1, t20, -t37, t21, t19, -t36, t25, t31 * t45 + t33 * t43, -t31 * t43 + t33 * t45 (-t28 - t29) * t17 * qJD(1), t24 + (t28 / 0.2e1 + t29 / 0.2e1) * t17 ^ 2, t20, t21, t37, t25, t36, t19, -t9 * qJD(3) + t10 * t44 (-t11 * t31 + t33 * t9) * qJD(1), t11 * qJD(3) - t10 * t42, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t19, -t37, -t36, t20, t21, t25, -t8 * qJD(3) + t7 * t42, t5 * qJD(3) + t7 * t44 (-t31 * t8 - t33 * t5) * qJD(1), t5 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t18, t12 ^ 2 / 0.2e1, -t12 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t6 * t12, t6 * t14 - t2 * t18, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg  = t13;
