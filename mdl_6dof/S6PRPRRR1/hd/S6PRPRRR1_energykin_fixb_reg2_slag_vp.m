% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:04
% EndTime: 2019-03-08 20:25:04
% DurationCPUTime: 0.16s
% Computational Cost: add. (418->53), mult. (956->134), div. (0->0), fcn. (689->12), ass. (0->47)
t49 = qJD(2) ^ 2;
t38 = t49 / 0.2e1;
t48 = cos(qJ(2));
t40 = sin(pkin(6));
t56 = qJD(1) * t40;
t28 = qJD(2) * pkin(2) + t48 * t56;
t39 = sin(pkin(12));
t41 = cos(pkin(12));
t46 = sin(qJ(2));
t52 = t46 * t56;
t19 = t39 * t28 + t41 * t52;
t17 = qJD(2) * pkin(8) + t19;
t42 = cos(pkin(6));
t34 = t42 * qJD(1) + qJD(3);
t47 = cos(qJ(4));
t32 = t47 * t34;
t45 = sin(qJ(4));
t10 = qJD(4) * pkin(4) + t32 + (-pkin(9) * qJD(2) - t17) * t45;
t13 = t47 * t17 + t45 * t34;
t54 = qJD(2) * t47;
t11 = pkin(9) * t54 + t13;
t44 = sin(qJ(5));
t58 = cos(qJ(5));
t6 = t44 * t10 + t58 * t11;
t57 = cos(qJ(6));
t55 = qJD(2) * t45;
t53 = qJD(2) * qJD(4);
t51 = qJD(2) * t56;
t18 = t41 * t28 - t39 * t52;
t24 = t44 * t55 - t58 * t54;
t5 = t58 * t10 - t44 * t11;
t14 = (-pkin(4) * t47 - pkin(3)) * qJD(2) - t18;
t50 = qJD(1) ^ 2;
t43 = sin(qJ(6));
t37 = qJD(4) + qJD(5);
t26 = (t44 * t47 + t58 * t45) * qJD(2);
t23 = qJD(6) + t24;
t22 = t57 * t26 + t43 * t37;
t20 = t43 * t26 - t57 * t37;
t16 = -qJD(2) * pkin(3) - t18;
t12 = -t45 * t17 + t32;
t7 = t24 * pkin(5) - t26 * pkin(10) + t14;
t4 = t37 * pkin(10) + t6;
t3 = -t37 * pkin(5) - t5;
t2 = t57 * t4 + t43 * t7;
t1 = -t43 * t4 + t57 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t50 / 0.2e1, 0, 0, 0, 0, 0, t38, t48 * t51, -t46 * t51, 0 (t42 ^ 2 / 0.2e1 + (t46 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1) * t40 ^ 2) * t50, 0, 0, 0, 0, 0, t38, t18 * qJD(2), -t19 * qJD(2), 0, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t45 ^ 2 * t38, t45 * t49 * t47, t45 * t53, t47 ^ 2 * t38, t47 * t53, qJD(4) ^ 2 / 0.2e1, t12 * qJD(4) - t16 * t54, -t13 * qJD(4) + t16 * t55 (-t12 * t45 + t13 * t47) * qJD(2), t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, t26 * t37, t24 ^ 2 / 0.2e1, -t24 * t37, t37 ^ 2 / 0.2e1, t14 * t24 + t5 * t37, t14 * t26 - t6 * t37, -t6 * t24 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t23, t20 ^ 2 / 0.2e1, -t20 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t3 * t20, -t2 * t23 + t3 * t22, -t1 * t22 - t2 * t20, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
