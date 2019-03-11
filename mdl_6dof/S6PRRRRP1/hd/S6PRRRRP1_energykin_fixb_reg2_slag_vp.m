% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:16
% EndTime: 2019-03-08 23:59:16
% DurationCPUTime: 0.15s
% Computational Cost: add. (499->58), mult. (1138->134), div. (0->0), fcn. (822->10), ass. (0->51)
t52 = qJD(2) ^ 2;
t63 = t52 / 0.2e1;
t50 = cos(qJ(3));
t51 = cos(qJ(2));
t44 = sin(pkin(6));
t60 = qJD(1) * t44;
t55 = t51 * t60;
t30 = -t55 + (-pkin(3) * t50 - pkin(2)) * qJD(2);
t47 = sin(qJ(4));
t57 = qJD(2) * t50;
t48 = sin(qJ(3));
t58 = qJD(2) * t48;
t62 = cos(qJ(4));
t32 = t47 * t58 - t62 * t57;
t34 = (t47 * t50 + t62 * t48) * qJD(2);
t14 = t32 * pkin(4) - t34 * pkin(10) + t30;
t46 = sin(qJ(5));
t61 = cos(qJ(5));
t49 = sin(qJ(2));
t36 = qJD(2) * pkin(8) + t49 * t60;
t45 = cos(pkin(6));
t59 = qJD(1) * t45;
t39 = t50 * t59;
t20 = qJD(3) * pkin(3) + t39 + (-pkin(9) * qJD(2) - t36) * t48;
t28 = t50 * t36 + t48 * t59;
t21 = pkin(9) * t57 + t28;
t10 = t47 * t20 + t62 * t21;
t43 = qJD(3) + qJD(4);
t8 = t43 * pkin(10) + t10;
t4 = t46 * t14 + t61 * t8;
t56 = qJD(2) * qJD(3);
t3 = t61 * t14 - t46 * t8;
t54 = qJD(2) * t60;
t9 = t62 * t20 - t47 * t21;
t7 = -t43 * pkin(4) - t9;
t53 = qJD(1) ^ 2;
t37 = -qJD(2) * pkin(2) - t55;
t31 = qJD(5) + t32;
t29 = t31 ^ 2 / 0.2e1;
t27 = -t48 * t36 + t39;
t26 = t61 * t34 + t46 * t43;
t24 = t46 * t34 - t61 * t43;
t23 = t26 ^ 2 / 0.2e1;
t22 = t24 ^ 2 / 0.2e1;
t16 = t26 * t31;
t15 = t24 * t31;
t13 = t26 * t24;
t5 = t24 * pkin(5) + qJD(6) + t7;
t2 = -t24 * qJ(6) + t4;
t1 = t31 * pkin(5) - t26 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t53 / 0.2e1, 0, 0, 0, 0, 0, t63, t51 * t54, -t49 * t54, 0 (t45 ^ 2 / 0.2e1 + (t49 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1) * t44 ^ 2) * t53, t48 ^ 2 * t63, t48 * t52 * t50, t48 * t56, t50 ^ 2 * t63, t50 * t56, qJD(3) ^ 2 / 0.2e1, t27 * qJD(3) - t37 * t57, -t28 * qJD(3) + t37 * t58 (-t27 * t48 + t28 * t50) * qJD(2), t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t43, t32 ^ 2 / 0.2e1, -t32 * t43, t43 ^ 2 / 0.2e1, t30 * t32 + t9 * t43, -t10 * t43 + t30 * t34, -t10 * t32 - t9 * t34, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t23, -t13, t16, t22, -t15, t29, t7 * t24 + t3 * t31, t7 * t26 - t4 * t31, -t4 * t24 - t3 * t26, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t23, -t13, t16, t22, -t15, t29, t1 * t31 + t5 * t24, -t2 * t31 + t5 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
