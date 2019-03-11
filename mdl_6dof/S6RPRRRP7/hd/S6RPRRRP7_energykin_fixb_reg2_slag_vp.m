% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP7
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:20:59
% EndTime: 2019-03-09 06:20:59
% DurationCPUTime: 0.18s
% Computational Cost: add. (788->61), mult. (1928->135), div. (0->0), fcn. (1413->8), ass. (0->49)
t53 = qJD(1) ^ 2;
t63 = t53 / 0.2e1;
t49 = sin(qJ(5));
t61 = cos(qJ(5));
t51 = sin(qJ(3));
t52 = cos(qJ(3));
t48 = cos(pkin(10));
t56 = qJD(1) * t48;
t47 = sin(pkin(10));
t57 = qJD(1) * t47;
t36 = t51 * t57 - t52 * t56;
t38 = (t47 * t52 + t48 * t51) * qJD(1);
t41 = qJD(2) + (-pkin(2) * t48 - pkin(1)) * qJD(1);
t20 = t36 * pkin(3) - t38 * pkin(8) + t41;
t58 = pkin(7) + qJ(2);
t39 = t58 * t57;
t40 = t58 * t56;
t25 = -t51 * t39 + t52 * t40;
t23 = qJD(3) * pkin(8) + t25;
t50 = sin(qJ(4));
t62 = cos(qJ(4));
t10 = t62 * t20 - t50 * t23;
t29 = t50 * qJD(3) + t62 * t38;
t32 = qJD(4) + t36;
t7 = t32 * pkin(4) - t29 * pkin(9) + t10;
t11 = t50 * t20 + t62 * t23;
t27 = -t62 * qJD(3) + t50 * t38;
t9 = -t27 * pkin(9) + t11;
t4 = t49 * t7 + t61 * t9;
t15 = t61 * t27 + t49 * t29;
t17 = -t49 * t27 + t61 * t29;
t60 = t17 * t15;
t31 = qJD(5) + t32;
t59 = t31 * t15;
t55 = t15 ^ 2 / 0.2e1;
t24 = -t52 * t39 - t51 * t40;
t3 = -t49 * t9 + t61 * t7;
t22 = -qJD(3) * pkin(3) - t24;
t13 = t27 * pkin(4) + t22;
t46 = t48 ^ 2;
t45 = t47 ^ 2;
t43 = -qJD(1) * pkin(1) + qJD(2);
t30 = t31 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t31;
t5 = t15 * pkin(5) - t17 * qJ(6) + t13;
t2 = t31 * qJ(6) + t4;
t1 = -t31 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t63, 0, 0, 0, 0, t45 * t63, t47 * t53 * t48, 0, t46 * t63, 0, 0, -t43 * t56, t43 * t57 (t45 + t46) * t53 * qJ(2), t43 ^ 2 / 0.2e1 + (t46 / 0.2e1 + t45 / 0.2e1) * qJ(2) ^ 2 * t53, t38 ^ 2 / 0.2e1, -t38 * t36, t38 * qJD(3), t36 ^ 2 / 0.2e1, -t36 * qJD(3), qJD(3) ^ 2 / 0.2e1, t24 * qJD(3) + t41 * t36, -t25 * qJD(3) + t41 * t38, -t24 * t38 - t25 * t36, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t32, t27 ^ 2 / 0.2e1, -t27 * t32, t32 ^ 2 / 0.2e1, t10 * t32 + t22 * t27, -t11 * t32 + t22 * t29, -t10 * t29 - t11 * t27, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t14, -t60, t12, t55, -t59, t30, t13 * t15 + t3 * t31, t13 * t17 - t4 * t31, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t14, t12, t60, t30, t59, t55, -t1 * t31 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t31, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
