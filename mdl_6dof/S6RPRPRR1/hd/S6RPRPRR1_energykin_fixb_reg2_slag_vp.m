% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:18
% EndTime: 2019-03-09 03:35:18
% DurationCPUTime: 0.16s
% Computational Cost: add. (657->61), mult. (1534->147), div. (0->0), fcn. (1047->10), ass. (0->48)
t51 = qJD(1) ^ 2;
t43 = t51 / 0.2e1;
t45 = sin(pkin(10));
t35 = (pkin(1) * t45 + pkin(7)) * qJD(1);
t50 = cos(qJ(3));
t40 = t50 * qJD(2);
t49 = sin(qJ(3));
t55 = qJ(4) * qJD(1);
t25 = qJD(3) * pkin(3) + t40 + (-t35 - t55) * t49;
t29 = t49 * qJD(2) + t50 * t35;
t26 = t50 * t55 + t29;
t44 = sin(pkin(11));
t58 = cos(pkin(11));
t12 = t58 * t25 - t44 * t26;
t33 = (t44 * t50 + t58 * t49) * qJD(1);
t10 = qJD(3) * pkin(4) - t33 * pkin(8) + t12;
t13 = t44 * t25 + t58 * t26;
t56 = qJD(1) * t50;
t57 = qJD(1) * t49;
t31 = t44 * t57 - t58 * t56;
t11 = -t31 * pkin(8) + t13;
t48 = sin(qJ(5));
t60 = cos(qJ(5));
t6 = t48 * t10 + t60 * t11;
t61 = pkin(1) * t51;
t59 = cos(qJ(6));
t54 = qJD(1) * qJD(3);
t46 = cos(pkin(10));
t53 = -pkin(1) * t46 - pkin(2);
t18 = t60 * t31 + t48 * t33;
t5 = t60 * t10 - t48 * t11;
t30 = qJD(4) + (-pkin(3) * t50 + t53) * qJD(1);
t21 = t31 * pkin(4) + t30;
t47 = sin(qJ(6));
t42 = qJD(3) ^ 2 / 0.2e1;
t41 = qJD(3) + qJD(5);
t36 = t53 * qJD(1);
t28 = -t49 * t35 + t40;
t20 = -t48 * t31 + t60 * t33;
t17 = qJD(6) + t18;
t16 = t59 * t20 + t47 * t41;
t14 = t47 * t20 - t59 * t41;
t7 = t18 * pkin(5) - t20 * pkin(9) + t21;
t4 = t41 * pkin(9) + t6;
t3 = -t41 * pkin(5) - t5;
t2 = t59 * t4 + t47 * t7;
t1 = -t47 * t4 + t59 * t7;
t8 = [0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t46 * t61, -t45 * t61, 0, qJD(2) ^ 2 / 0.2e1 + (t45 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t51, t49 ^ 2 * t43, t49 * t51 * t50, t49 * t54, t50 ^ 2 * t43, t50 * t54, t42, t28 * qJD(3) - t36 * t56, -t29 * qJD(3) + t36 * t57 (-t28 * t49 + t29 * t50) * qJD(1), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * qJD(3), t31 ^ 2 / 0.2e1, -t31 * qJD(3), t42, t12 * qJD(3) + t30 * t31, -t13 * qJD(3) + t30 * t33, -t12 * t33 - t13 * t31, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t41, t18 ^ 2 / 0.2e1, -t18 * t41, t41 ^ 2 / 0.2e1, t21 * t18 + t5 * t41, t21 * t20 - t6 * t41, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t17, t14 ^ 2 / 0.2e1, -t14 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t14, t3 * t16 - t2 * t17, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
