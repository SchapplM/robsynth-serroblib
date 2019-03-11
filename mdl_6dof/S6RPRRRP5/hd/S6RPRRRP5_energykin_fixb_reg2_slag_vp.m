% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:24
% EndTime: 2019-03-09 06:12:24
% DurationCPUTime: 0.19s
% Computational Cost: add. (797->61), mult. (2037->135), div. (0->0), fcn. (1528->8), ass. (0->49)
t52 = qJD(1) ^ 2;
t63 = t52 / 0.2e1;
t51 = sin(qJ(3));
t48 = cos(pkin(10));
t55 = qJD(1) * t48;
t47 = sin(pkin(10));
t56 = qJD(1) * t47;
t62 = cos(qJ(3));
t35 = t51 * t56 - t62 * t55;
t37 = (t62 * t47 + t48 * t51) * qJD(1);
t50 = sin(qJ(4));
t61 = cos(qJ(4));
t25 = t61 * t35 + t50 * t37;
t27 = -t50 * t35 + t61 * t37;
t40 = qJD(2) + (-pkin(2) * t48 - pkin(1)) * qJD(1);
t30 = t35 * pkin(3) + t40;
t12 = t25 * pkin(4) - t27 * pkin(9) + t30;
t49 = sin(qJ(5));
t60 = cos(qJ(5));
t57 = pkin(7) + qJ(2);
t38 = t57 * t56;
t39 = t57 * t55;
t28 = -t62 * t38 - t51 * t39;
t17 = qJD(3) * pkin(3) - t37 * pkin(8) + t28;
t29 = -t51 * t38 + t62 * t39;
t19 = -t35 * pkin(8) + t29;
t11 = t50 * t17 + t61 * t19;
t46 = qJD(3) + qJD(4);
t8 = t46 * pkin(9) + t11;
t4 = t49 * t12 + t60 * t8;
t20 = t49 * t27 - t60 * t46;
t22 = t60 * t27 + t49 * t46;
t59 = t22 * t20;
t24 = qJD(5) + t25;
t58 = t24 * t20;
t54 = t20 ^ 2 / 0.2e1;
t10 = t61 * t17 - t50 * t19;
t3 = t60 * t12 - t49 * t8;
t7 = -t46 * pkin(4) - t10;
t45 = t48 ^ 2;
t44 = t47 ^ 2;
t43 = -qJD(1) * pkin(1) + qJD(2);
t23 = t24 ^ 2 / 0.2e1;
t18 = t22 ^ 2 / 0.2e1;
t13 = t22 * t24;
t5 = t20 * pkin(5) - t22 * qJ(6) + t7;
t2 = t24 * qJ(6) + t4;
t1 = -t24 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t63, 0, 0, 0, 0, t44 * t63, t47 * t52 * t48, 0, t45 * t63, 0, 0, -t43 * t55, t43 * t56 (t44 + t45) * t52 * qJ(2), t43 ^ 2 / 0.2e1 + (t45 / 0.2e1 + t44 / 0.2e1) * qJ(2) ^ 2 * t52, t37 ^ 2 / 0.2e1, -t37 * t35, t37 * qJD(3), t35 ^ 2 / 0.2e1, -t35 * qJD(3), qJD(3) ^ 2 / 0.2e1, t28 * qJD(3) + t40 * t35, -t29 * qJD(3) + t40 * t37, -t28 * t37 - t29 * t35, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t46, t25 ^ 2 / 0.2e1, -t25 * t46, t46 ^ 2 / 0.2e1, t10 * t46 + t30 * t25, -t11 * t46 + t30 * t27, -t10 * t27 - t11 * t25, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t18, -t59, t13, t54, -t58, t23, t7 * t20 + t3 * t24, t7 * t22 - t4 * t24, -t4 * t20 - t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t18, t13, t59, t23, t58, t54, -t1 * t24 + t5 * t20, t1 * t22 - t2 * t20, t2 * t24 - t5 * t22, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
