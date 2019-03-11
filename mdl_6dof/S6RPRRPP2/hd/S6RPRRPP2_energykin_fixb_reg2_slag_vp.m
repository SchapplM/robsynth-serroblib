% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:20
% EndTime: 2019-03-09 04:33:20
% DurationCPUTime: 0.15s
% Computational Cost: add. (317->51), mult. (702->111), div. (0->0), fcn. (389->6), ass. (0->41)
t41 = qJD(1) ^ 2;
t35 = t41 / 0.2e1;
t53 = pkin(4) + pkin(5);
t52 = pkin(1) * t41;
t51 = cos(qJ(4));
t38 = sin(qJ(4));
t39 = sin(qJ(3));
t49 = qJD(1) * t39;
t23 = -t51 * qJD(3) + t38 * t49;
t40 = cos(qJ(3));
t48 = t40 * qJD(1);
t30 = -qJD(4) + t48;
t50 = t23 * t30;
t25 = t38 * qJD(3) + t51 * t49;
t9 = t25 * t23;
t16 = t25 * t30;
t36 = sin(pkin(9));
t26 = (pkin(1) * t36 + pkin(7)) * qJD(1);
t18 = t39 * qJD(2) + t40 * t26;
t14 = qJD(3) * pkin(8) + t18;
t37 = cos(pkin(9));
t46 = -pkin(1) * t37 - pkin(2);
t15 = (-pkin(3) * t40 - pkin(8) * t39 + t46) * qJD(1);
t8 = t51 * t14 + t38 * t15;
t17 = t40 * qJD(2) - t39 * t26;
t19 = t23 ^ 2 / 0.2e1;
t28 = t30 ^ 2 / 0.2e1;
t47 = qJD(1) * qJD(3);
t5 = -t30 * qJ(5) + t8;
t45 = qJD(3) * pkin(3) + t17;
t7 = -t38 * t14 + t51 * t15;
t44 = qJD(5) - t7;
t43 = t25 * qJ(5) + t45;
t27 = t46 * qJD(1);
t20 = t25 ^ 2 / 0.2e1;
t6 = t23 * pkin(4) - t43;
t4 = t30 * pkin(4) + t44;
t3 = -t53 * t23 + qJD(6) + t43;
t2 = t23 * qJ(6) + t5;
t1 = -t25 * qJ(6) + t53 * t30 + t44;
t10 = [0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t37 * t52, -t36 * t52, 0, qJD(2) ^ 2 / 0.2e1 + (t36 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t41, t39 ^ 2 * t35, t39 * t41 * t40, t39 * t47, t40 ^ 2 * t35, t40 * t47, qJD(3) ^ 2 / 0.2e1, t17 * qJD(3) - t27 * t48, -t18 * qJD(3) + t27 * t49 (-t17 * t39 + t18 * t40) * qJD(1), t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t20, -t9, -t16, t19, t50, t28, -t23 * t45 - t7 * t30, -t25 * t45 + t8 * t30, -t8 * t23 - t7 * t25, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1, t20, -t16, t9, t28, -t50, t19, t6 * t23 + t4 * t30, -t5 * t23 + t4 * t25, -t6 * t25 - t5 * t30, t5 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t20, t9, t16, t19, t50, t28, t1 * t30 - t3 * t23, -t2 * t30 + t3 * t25, -t1 * t25 + t2 * t23, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t10;
