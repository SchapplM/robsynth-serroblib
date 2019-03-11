% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:45
% EndTime: 2019-03-08 20:48:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (351->52), mult. (734->122), div. (0->0), fcn. (475->10), ass. (0->47)
t44 = qJD(2) ^ 2;
t35 = t44 / 0.2e1;
t45 = qJD(1) ^ 2;
t58 = t45 / 0.2e1;
t43 = cos(qJ(2));
t36 = sin(pkin(6));
t55 = qJD(1) * t36;
t46 = -t43 * t55 + qJD(3);
t21 = (-pkin(2) - pkin(8)) * qJD(2) + t46;
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t37 = cos(pkin(6));
t54 = qJD(1) * t37;
t15 = t40 * t21 + t42 * t54;
t13 = qJD(4) * pkin(9) + t15;
t41 = sin(qJ(2));
t48 = t41 * t55;
t18 = t48 + (pkin(4) * t40 - pkin(9) * t42 + qJ(3)) * qJD(2);
t39 = sin(qJ(5));
t57 = cos(qJ(5));
t6 = t57 * t13 + t39 * t18;
t56 = cos(qJ(6));
t53 = qJD(2) * t42;
t26 = qJD(2) * qJ(3) + t48;
t52 = t26 * qJD(2);
t51 = t40 * qJD(2);
t50 = t26 ^ 2 / 0.2e1;
t49 = qJD(2) * qJD(4);
t47 = qJD(2) * t55;
t5 = -t39 * t13 + t57 * t18;
t14 = t42 * t21 - t40 * t54;
t31 = qJD(5) + t51;
t12 = -qJD(4) * pkin(4) - t14;
t38 = sin(qJ(6));
t32 = t37 ^ 2 * t58;
t28 = qJD(6) + t31;
t25 = t39 * qJD(4) + t57 * t53;
t23 = -t57 * qJD(4) + t39 * t53;
t22 = -qJD(2) * pkin(2) + t46;
t11 = -t38 * t23 + t56 * t25;
t9 = t56 * t23 + t38 * t25;
t7 = t23 * pkin(5) + t12;
t4 = -t23 * pkin(10) + t6;
t3 = t31 * pkin(5) - t25 * pkin(10) + t5;
t2 = t38 * t3 + t56 * t4;
t1 = t56 * t3 - t38 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, t35, t43 * t47, -t41 * t47, 0, t32 + (t41 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1) * t45 * t36 ^ 2, t35, 0, 0, 0, 0, 0, 0, t22 * qJD(2), t52, t32 + t50 + t22 ^ 2 / 0.2e1, t42 ^ 2 * t35, -t42 * t44 * t40, t42 * t49, t40 ^ 2 * t35, -t40 * t49, qJD(4) ^ 2 / 0.2e1, t14 * qJD(4) + t26 * t51, -t15 * qJD(4) + t42 * t52 (-t14 * t42 - t15 * t40) * qJD(2), t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t50, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t31, t23 ^ 2 / 0.2e1, -t23 * t31, t31 ^ 2 / 0.2e1, t12 * t23 + t5 * t31, t12 * t25 - t6 * t31, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t9, t11 * t28, t9 ^ 2 / 0.2e1, -t9 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t9, t7 * t11 - t2 * t28, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
