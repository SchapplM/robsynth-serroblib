% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:10:56
% EndTime: 2019-03-08 21:10:56
% DurationCPUTime: 0.14s
% Computational Cost: add. (235->52), mult. (559->117), div. (0->0), fcn. (311->8), ass. (0->46)
t44 = qJD(2) ^ 2;
t56 = t44 / 0.2e1;
t55 = -pkin(3) - pkin(4);
t40 = sin(qJ(2));
t35 = sin(pkin(6));
t54 = qJD(1) * t35;
t19 = qJD(2) * pkin(8) + t40 * t54;
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t36 = cos(pkin(6));
t53 = qJD(1) * t36;
t12 = t42 * t19 + t39 * t53;
t43 = cos(qJ(2));
t20 = -qJD(2) * pkin(2) - t43 * t54;
t52 = qJD(2) * t42;
t51 = t39 * qJD(2);
t50 = qJ(5) * qJD(2);
t49 = qJD(2) * qJD(3);
t10 = qJD(3) * qJ(4) + t12;
t48 = qJD(2) * t54;
t13 = -pkin(3) * t52 - qJ(4) * t51 + t20;
t11 = -t39 * t19 + t42 * t53;
t8 = pkin(4) * t52 + qJD(5) - t13;
t47 = qJD(4) - t11;
t7 = t42 * t50 - t10;
t46 = -t39 * t50 + t47;
t45 = qJD(1) ^ 2;
t41 = cos(qJ(6));
t38 = sin(qJ(6));
t33 = qJD(3) ^ 2 / 0.2e1;
t28 = t42 * t49;
t27 = t39 * t49;
t26 = t42 ^ 2 * t56;
t25 = t39 ^ 2 * t56;
t24 = qJD(6) + t51;
t23 = t39 * t44 * t42;
t17 = t38 * qJD(3) + t41 * t52;
t16 = -t41 * qJD(3) + t38 * t52;
t9 = -qJD(3) * pkin(3) + t47;
t6 = qJD(3) * pkin(5) - t7;
t5 = t55 * qJD(3) + t46;
t4 = (-pkin(9) + t55) * qJD(3) + t46;
t3 = (pkin(5) * t39 + pkin(9) * t42) * qJD(2) + t8;
t2 = t38 * t3 + t41 * t4;
t1 = t41 * t3 - t38 * t4;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t45 / 0.2e1, 0, 0, 0, 0, 0, t56, t43 * t48, -t40 * t48, 0 (t36 ^ 2 / 0.2e1 + (t40 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1) * t35 ^ 2) * t45, t25, t23, t27, t26, t28, t33, t11 * qJD(3) - t20 * t52, -t12 * qJD(3) + t20 * t51 (-t11 * t39 + t12 * t42) * qJD(2), t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t25, t27, -t23, t33, -t28, t26, -t9 * qJD(3) - t13 * t52 (t10 * t42 + t39 * t9) * qJD(2), t10 * qJD(3) - t13 * t51, t10 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t26, t23, t28, t25, t27, t33, -t7 * qJD(3) + t8 * t51, t5 * qJD(3) - t8 * t52 (-t39 * t5 + t42 * t7) * qJD(2), t5 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t16, -t17 * t24, t16 ^ 2 / 0.2e1, t16 * t24, t24 ^ 2 / 0.2e1, t1 * t24 - t6 * t16, -t6 * t17 - t2 * t24, t1 * t17 + t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg  = t14;
