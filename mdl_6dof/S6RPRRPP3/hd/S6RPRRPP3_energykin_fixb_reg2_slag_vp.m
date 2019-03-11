% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPP3
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:36:57
% EndTime: 2019-03-09 04:36:57
% DurationCPUTime: 0.16s
% Computational Cost: add. (317->53), mult. (702->111), div. (0->0), fcn. (389->6), ass. (0->41)
t40 = qJD(1) ^ 2;
t34 = t40 / 0.2e1;
t51 = pkin(1) * t40;
t50 = cos(qJ(4));
t37 = sin(qJ(4));
t38 = sin(qJ(3));
t47 = qJD(1) * t38;
t23 = -t50 * qJD(3) + t37 * t47;
t25 = t37 * qJD(3) + t50 * t47;
t49 = t23 * t25;
t39 = cos(qJ(3));
t46 = t39 * qJD(1);
t30 = -qJD(4) + t46;
t15 = t25 * t30;
t16 = t30 * t23;
t48 = pkin(4) + qJ(6);
t35 = sin(pkin(9));
t27 = (pkin(1) * t35 + pkin(7)) * qJD(1);
t18 = t38 * qJD(2) + t39 * t27;
t13 = qJD(3) * pkin(8) + t18;
t36 = cos(pkin(9));
t44 = -pkin(1) * t36 - pkin(2);
t14 = (-pkin(3) * t39 - pkin(8) * t38 + t44) * qJD(1);
t8 = t50 * t13 + t37 * t14;
t19 = t23 ^ 2 / 0.2e1;
t20 = t25 ^ 2 / 0.2e1;
t45 = qJD(1) * qJD(3);
t17 = t39 * qJD(2) - t38 * t27;
t5 = t30 * qJ(5) - t8;
t7 = -t37 * t13 + t50 * t14;
t43 = qJD(5) - t7;
t12 = -qJD(3) * pkin(3) - t17;
t42 = -t25 * qJ(5) + t12;
t29 = t30 ^ 2 / 0.2e1;
t28 = t44 * qJD(1);
t6 = t23 * pkin(4) + t42;
t4 = t30 * pkin(4) + t43;
t3 = t48 * t23 + t42;
t2 = -t23 * pkin(5) + qJD(6) - t5;
t1 = t25 * pkin(5) + t48 * t30 + t43;
t9 = [0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t36 * t51, -t35 * t51, 0, qJD(2) ^ 2 / 0.2e1 + (t35 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t40, t38 ^ 2 * t34, t38 * t40 * t39, t38 * t45, t39 ^ 2 * t34, t39 * t45, qJD(3) ^ 2 / 0.2e1, t17 * qJD(3) - t28 * t46, -t18 * qJD(3) + t28 * t47 (-t17 * t38 + t18 * t39) * qJD(1), t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t20, -t49, -t15, t19, t16, t29, t12 * t23 - t7 * t30, t12 * t25 + t8 * t30, -t8 * t23 - t7 * t25, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t29, t15, -t16, t20, -t49, t19, t5 * t23 + t4 * t25, -t6 * t23 - t4 * t30, -t6 * t25 + t5 * t30, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t29, -t16, -t15, t19, t49, t20, t1 * t25 - t2 * t23, -t2 * t30 - t3 * t25, t1 * t30 + t3 * t23, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg  = t9;
