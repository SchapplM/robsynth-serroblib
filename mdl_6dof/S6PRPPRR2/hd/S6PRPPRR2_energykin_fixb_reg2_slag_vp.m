% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:00
% EndTime: 2019-03-08 19:20:00
% DurationCPUTime: 0.13s
% Computational Cost: add. (219->43), mult. (494->102), div. (0->0), fcn. (317->10), ass. (0->41)
t36 = qJD(2) ^ 2;
t26 = t36 / 0.2e1;
t30 = cos(pkin(6));
t22 = t30 * qJD(1) + qJD(3);
t32 = sin(qJ(5));
t34 = cos(qJ(5));
t35 = cos(qJ(2));
t28 = sin(pkin(6));
t46 = qJD(1) * t28;
t18 = qJD(2) * pkin(2) + t35 * t46;
t27 = sin(pkin(11));
t29 = cos(pkin(11));
t33 = sin(qJ(2));
t40 = t33 * t46;
t13 = t29 * t18 - t27 * t40;
t38 = qJD(4) - t13;
t9 = (-pkin(3) - pkin(8)) * qJD(2) + t38;
t6 = t34 * t22 + t32 * t9;
t47 = cos(qJ(6));
t45 = qJD(2) * t34;
t14 = t27 * t18 + t29 * t40;
t11 = qJD(2) * qJ(4) + t14;
t44 = t11 * qJD(2);
t43 = t32 * qJD(2);
t42 = t11 ^ 2 / 0.2e1;
t41 = qJD(2) * qJD(5);
t39 = qJD(2) * t46;
t5 = -t32 * t22 + t34 * t9;
t37 = qJD(1) ^ 2;
t31 = sin(qJ(6));
t23 = qJD(6) + t43;
t21 = t22 ^ 2 / 0.2e1;
t17 = t31 * qJD(5) + t47 * t45;
t15 = -t47 * qJD(5) + t31 * t45;
t10 = -qJD(2) * pkin(3) + t38;
t7 = (pkin(5) * t32 - pkin(9) * t34 + qJ(4)) * qJD(2) + t14;
t4 = qJD(5) * pkin(9) + t6;
t3 = -qJD(5) * pkin(5) - t5;
t2 = t31 * t7 + t47 * t4;
t1 = -t31 * t4 + t47 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t37 / 0.2e1, 0, 0, 0, 0, 0, t26, t35 * t39, -t33 * t39, 0 (t30 ^ 2 / 0.2e1 + (t33 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1) * t28 ^ 2) * t37, 0, 0, 0, 0, 0, t26, t13 * qJD(2), -t14 * qJD(2), 0, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t21, t26, 0, 0, 0, 0, 0, 0, t10 * qJD(2), t44, t21 + t42 + t10 ^ 2 / 0.2e1, t34 ^ 2 * t26, -t34 * t36 * t32, t34 * t41, t32 ^ 2 * t26, -t32 * t41, qJD(5) ^ 2 / 0.2e1, t5 * qJD(5) + t11 * t43, -t6 * qJD(5) + t34 * t44 (-t32 * t6 - t34 * t5) * qJD(2), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t42, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t23, t15 ^ 2 / 0.2e1, -t15 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t3 * t15, t3 * t17 - t2 * t23, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
