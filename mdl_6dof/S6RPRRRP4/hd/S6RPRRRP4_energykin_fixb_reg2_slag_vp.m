% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP4
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:44
% EndTime: 2019-03-09 06:08:44
% DurationCPUTime: 0.19s
% Computational Cost: add. (800->63), mult. (2043->135), div. (0->0), fcn. (1534->8), ass. (0->49)
t56 = qJD(1) ^ 2;
t64 = t56 / 0.2e1;
t55 = sin(qJ(3));
t52 = cos(pkin(10));
t58 = qJD(1) * t52;
t51 = sin(pkin(10));
t59 = qJD(1) * t51;
t63 = cos(qJ(3));
t39 = t55 * t59 - t63 * t58;
t41 = (t63 * t51 + t52 * t55) * qJD(1);
t54 = sin(qJ(4));
t62 = cos(qJ(4));
t29 = t62 * t39 + t54 * t41;
t31 = -t54 * t39 + t62 * t41;
t44 = qJD(2) + (-pkin(2) * t52 - pkin(1)) * qJD(1);
t34 = t39 * pkin(3) + t44;
t13 = t29 * pkin(4) - t31 * pkin(9) + t34;
t53 = sin(qJ(5));
t61 = cos(qJ(5));
t60 = pkin(7) + qJ(2);
t42 = t60 * t59;
t43 = t60 * t58;
t32 = -t63 * t42 - t55 * t43;
t20 = qJD(3) * pkin(3) - t41 * pkin(8) + t32;
t33 = -t55 * t42 + t63 * t43;
t23 = -t39 * pkin(8) + t33;
t12 = t54 * t20 + t62 * t23;
t50 = qJD(3) + qJD(4);
t8 = t50 * pkin(9) + t12;
t4 = t53 * t13 + t61 * t8;
t3 = t61 * t13 - t53 * t8;
t11 = t62 * t20 - t54 * t23;
t7 = -t50 * pkin(4) - t11;
t49 = t52 ^ 2;
t48 = t51 ^ 2;
t47 = -qJD(1) * pkin(1) + qJD(2);
t28 = qJD(5) + t29;
t27 = t28 ^ 2 / 0.2e1;
t26 = t61 * t31 + t53 * t50;
t24 = t53 * t31 - t61 * t50;
t22 = t26 ^ 2 / 0.2e1;
t21 = t24 ^ 2 / 0.2e1;
t16 = t26 * t28;
t15 = t24 * t28;
t14 = t26 * t24;
t5 = t24 * pkin(5) + qJD(6) + t7;
t2 = -t24 * qJ(6) + t4;
t1 = t28 * pkin(5) - t26 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t64, 0, 0, 0, 0, t48 * t64, t51 * t56 * t52, 0, t49 * t64, 0, 0, -t47 * t58, t47 * t59 (t48 + t49) * t56 * qJ(2), t47 ^ 2 / 0.2e1 + (t49 / 0.2e1 + t48 / 0.2e1) * qJ(2) ^ 2 * t56, t41 ^ 2 / 0.2e1, -t39 * t41, qJD(3) * t41, t39 ^ 2 / 0.2e1, -t39 * qJD(3), qJD(3) ^ 2 / 0.2e1, t32 * qJD(3) + t44 * t39, -t33 * qJD(3) + t44 * t41, -t32 * t41 - t33 * t39, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t50, t29 ^ 2 / 0.2e1, -t29 * t50, t50 ^ 2 / 0.2e1, t11 * t50 + t34 * t29, -t12 * t50 + t34 * t31, -t11 * t31 - t12 * t29, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t22, -t14, t16, t21, -t15, t27, t7 * t24 + t3 * t28, t7 * t26 - t4 * t28, -t4 * t24 - t3 * t26, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t22, -t14, t16, t21, -t15, t27, t1 * t28 + t5 * t24, -t2 * t28 + t5 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
