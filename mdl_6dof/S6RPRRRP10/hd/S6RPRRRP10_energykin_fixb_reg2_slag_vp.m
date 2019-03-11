% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:31
% EndTime: 2019-03-09 06:32:31
% DurationCPUTime: 0.14s
% Computational Cost: add. (480->53), mult. (940->115), div. (0->0), fcn. (558->6), ass. (0->43)
t43 = qJD(1) ^ 2;
t36 = t43 / 0.2e1;
t39 = sin(qJ(5));
t51 = cos(qJ(5));
t41 = sin(qJ(3));
t42 = cos(qJ(3));
t22 = (pkin(3) * t41 - pkin(8) * t42 + qJ(2)) * qJD(1);
t31 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t23 = qJD(3) * pkin(8) + t41 * t31;
t40 = sin(qJ(4));
t52 = cos(qJ(4));
t11 = t52 * t22 - t40 * t23;
t47 = qJD(1) * t42;
t27 = t40 * qJD(3) + t52 * t47;
t32 = t41 * qJD(1) + qJD(4);
t7 = t32 * pkin(4) - t27 * pkin(9) + t11;
t12 = t40 * t22 + t52 * t23;
t25 = -t52 * qJD(3) + t40 * t47;
t9 = -t25 * pkin(9) + t12;
t4 = t39 * t7 + t51 * t9;
t14 = t51 * t25 + t39 * t27;
t16 = -t39 * t25 + t51 * t27;
t50 = t16 * t14;
t30 = qJD(5) + t32;
t49 = t30 * t14;
t48 = t43 * qJ(2);
t46 = qJD(3) * t31;
t45 = t14 ^ 2 / 0.2e1;
t44 = qJD(1) * qJD(3);
t24 = -qJD(3) * pkin(3) - t42 * t31;
t3 = -t39 * t9 + t51 * t7;
t17 = t25 * pkin(4) + t24;
t38 = t42 ^ 2;
t37 = t41 ^ 2;
t34 = qJ(2) ^ 2 * t36;
t33 = -pkin(1) * qJD(1) + qJD(2);
t28 = t30 ^ 2 / 0.2e1;
t13 = t16 ^ 2 / 0.2e1;
t10 = t16 * t30;
t5 = t14 * pkin(5) - t16 * qJ(6) + t17;
t2 = t30 * qJ(6) + t4;
t1 = -t30 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t36, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, t33 * qJD(1), t48, t34 + t33 ^ 2 / 0.2e1, t38 * t36, -t42 * t43 * t41, t42 * t44, t37 * t36, -t41 * t44, qJD(3) ^ 2 / 0.2e1, t41 * t48 + t42 * t46, -t41 * t46 + t42 * t48 (-t37 - t38) * t31 * qJD(1), t34 + (t37 / 0.2e1 + t38 / 0.2e1) * t31 ^ 2, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t32, t25 ^ 2 / 0.2e1, -t25 * t32, t32 ^ 2 / 0.2e1, t11 * t32 + t24 * t25, -t12 * t32 + t24 * t27, -t11 * t27 - t12 * t25, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t13, -t50, t10, t45, -t49, t28, t17 * t14 + t3 * t30, t17 * t16 - t4 * t30, -t4 * t14 - t3 * t16, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t13, t10, t50, t28, t49, t45, -t1 * t30 + t5 * t14, t1 * t16 - t2 * t14, -t5 * t16 + t2 * t30, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
