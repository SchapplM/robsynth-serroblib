% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:00
% EndTime: 2019-03-09 03:13:00
% DurationCPUTime: 0.16s
% Computational Cost: add. (283->51), mult. (622->109), div. (0->0), fcn. (310->6), ass. (0->47)
t42 = qJD(1) ^ 2;
t34 = t42 / 0.2e1;
t56 = -pkin(3) - pkin(8);
t41 = cos(qJ(3));
t39 = sin(qJ(3));
t37 = cos(pkin(9));
t48 = -pkin(1) * t37 - pkin(2);
t44 = -qJ(4) * t39 + t48;
t10 = (t56 * t41 + t44) * qJD(1);
t38 = sin(qJ(5));
t40 = cos(qJ(5));
t36 = sin(pkin(9));
t23 = (pkin(1) * t36 + pkin(7)) * qJD(1);
t15 = t41 * qJD(2) - t39 * t23;
t45 = qJD(4) - t15;
t51 = t39 * qJD(1);
t7 = pkin(4) * t51 + t56 * qJD(3) + t45;
t4 = t40 * t10 + t38 * t7;
t55 = pkin(1) * t42;
t52 = qJD(1) * t41;
t20 = t38 * qJD(3) + t40 * t52;
t22 = t40 * qJD(3) - t38 * t52;
t54 = t22 * t20;
t28 = qJD(5) + t51;
t53 = t28 * t20;
t16 = t39 * qJD(2) + t41 * t23;
t50 = t20 ^ 2 / 0.2e1;
t49 = qJD(1) * qJD(3);
t12 = -qJD(3) * qJ(4) - t16;
t47 = t39 * t49;
t46 = t41 * t49;
t9 = pkin(4) * t52 - t12;
t3 = -t38 * t10 + t40 * t7;
t33 = qJD(3) ^ 2 / 0.2e1;
t30 = t41 ^ 2 * t34;
t29 = t39 ^ 2 * t34;
t27 = t39 * t42 * t41;
t25 = t28 ^ 2 / 0.2e1;
t24 = t48 * qJD(1);
t17 = t22 ^ 2 / 0.2e1;
t14 = t22 * t28;
t13 = (-pkin(3) * t41 + t44) * qJD(1);
t11 = -qJD(3) * pkin(3) + t45;
t5 = t20 * pkin(5) - t22 * qJ(6) + t9;
t2 = t28 * qJ(6) + t4;
t1 = -t28 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t37 * t55, -t36 * t55, 0, qJD(2) ^ 2 / 0.2e1 + (t36 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t42, t29, t27, t47, t30, t46, t33, t15 * qJD(3) - t24 * t52, -t16 * qJD(3) + t24 * t51 (-t15 * t39 + t16 * t41) * qJD(1), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t33, -t47, -t46, t29, t27, t30 (t11 * t39 - t12 * t41) * qJD(1), t11 * qJD(3) + t13 * t52, -t12 * qJD(3) - t13 * t51, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t17, -t54, t14, t50, -t53, t25, t9 * t20 + t3 * t28, t9 * t22 - t4 * t28, -t4 * t20 - t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t17, t14, t54, t25, t53, t50, -t1 * t28 + t5 * t20, t1 * t22 - t2 * t20, t2 * t28 - t5 * t22, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
