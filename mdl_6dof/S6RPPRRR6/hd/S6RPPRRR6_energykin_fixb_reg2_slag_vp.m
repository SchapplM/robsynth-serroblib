% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:28
% EndTime: 2019-03-09 02:31:28
% DurationCPUTime: 0.14s
% Computational Cost: add. (311->47), mult. (574->110), div. (0->0), fcn. (293->6), ass. (0->38)
t38 = qJD(1) ^ 2;
t30 = t38 / 0.2e1;
t47 = cos(qJ(5));
t46 = cos(qJ(6));
t45 = -pkin(1) - qJ(3);
t36 = sin(qJ(4));
t37 = cos(qJ(4));
t13 = -qJD(2) + (pkin(4) * t36 - pkin(8) * t37 - t45) * qJD(1);
t27 = qJ(2) * qJD(1) + qJD(3);
t22 = -pkin(7) * qJD(1) + t27;
t16 = qJD(4) * pkin(8) + t36 * t22;
t35 = sin(qJ(5));
t6 = t35 * t13 + t47 * t16;
t44 = qJD(1) * t37;
t43 = qJD(4) * t22;
t23 = -t45 * qJD(1) - qJD(2);
t42 = t23 * qJD(1);
t41 = t36 * qJD(1);
t40 = t23 ^ 2 / 0.2e1;
t39 = qJD(1) * qJD(4);
t5 = t47 * t13 - t35 * t16;
t26 = qJD(5) + t41;
t17 = -qJD(4) * pkin(4) - t37 * t22;
t34 = sin(qJ(6));
t33 = t37 ^ 2;
t32 = t36 ^ 2;
t28 = -qJD(1) * pkin(1) + qJD(2);
t25 = qJD(6) + t26;
t20 = t35 * qJD(4) + t47 * t44;
t18 = -t47 * qJD(4) + t35 * t44;
t10 = t18 * pkin(5) + t17;
t9 = -t34 * t18 + t46 * t20;
t7 = t46 * t18 + t34 * t20;
t4 = -t18 * pkin(9) + t6;
t3 = t26 * pkin(5) - t20 * pkin(9) + t5;
t2 = t34 * t3 + t46 * t4;
t1 = t46 * t3 - t34 * t4;
t8 = [0, 0, 0, 0, 0, t30, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, t28 * qJD(1), t38 * qJ(2), qJ(2) ^ 2 * t30 + t28 ^ 2 / 0.2e1, t30, 0, 0, 0, 0, 0, 0, t27 * qJD(1), t42, t40 + t27 ^ 2 / 0.2e1, t33 * t30, -t37 * t38 * t36, t37 * t39, t32 * t30, -t36 * t39, qJD(4) ^ 2 / 0.2e1, t23 * t41 + t37 * t43, -t36 * t43 + t37 * t42 (-t32 - t33) * t22 * qJD(1), t40 + (t32 / 0.2e1 + t33 / 0.2e1) * t22 ^ 2, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t26, t18 ^ 2 / 0.2e1, -t18 * t26, t26 ^ 2 / 0.2e1, t17 * t18 + t5 * t26, t17 * t20 - t6 * t26, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t25, t7 ^ 2 / 0.2e1, -t7 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t10 * t7, t10 * t9 - t2 * t25, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t8;
