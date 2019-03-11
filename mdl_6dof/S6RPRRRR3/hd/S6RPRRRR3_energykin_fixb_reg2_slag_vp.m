% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:28
% EndTime: 2019-03-09 07:02:28
% DurationCPUTime: 0.17s
% Computational Cost: add. (721->60), mult. (1547->148), div. (0->0), fcn. (1021->10), ass. (0->46)
t55 = qJD(1) ^ 2;
t47 = t55 / 0.2e1;
t64 = pkin(1) * t55;
t48 = sin(pkin(11));
t39 = (pkin(1) * t48 + pkin(7)) * qJD(1);
t53 = sin(qJ(3));
t54 = cos(qJ(3));
t30 = t53 * qJD(2) + t54 * t39;
t27 = qJD(3) * pkin(8) + t30;
t49 = cos(pkin(11));
t57 = -pkin(1) * t49 - pkin(2);
t28 = (-pkin(3) * t54 - pkin(8) * t53 + t57) * qJD(1);
t52 = sin(qJ(4));
t63 = cos(qJ(4));
t16 = -t52 * t27 + t63 * t28;
t60 = qJD(1) * t53;
t36 = t52 * qJD(3) + t63 * t60;
t59 = t54 * qJD(1);
t43 = -qJD(4) + t59;
t12 = -t43 * pkin(4) - t36 * pkin(9) + t16;
t17 = t63 * t27 + t52 * t28;
t34 = -t63 * qJD(3) + t52 * t60;
t15 = -t34 * pkin(9) + t17;
t51 = sin(qJ(5));
t62 = cos(qJ(5));
t6 = t51 * t12 + t62 * t15;
t61 = cos(qJ(6));
t58 = qJD(1) * qJD(3);
t5 = t62 * t12 - t51 * t15;
t29 = t54 * qJD(2) - t53 * t39;
t41 = -qJD(5) + t43;
t26 = -qJD(3) * pkin(3) - t29;
t19 = t34 * pkin(4) + t26;
t50 = sin(qJ(6));
t40 = t57 * qJD(1);
t37 = -qJD(6) + t41;
t22 = -t51 * t34 + t62 * t36;
t20 = t62 * t34 + t51 * t36;
t13 = t20 * pkin(5) + t19;
t11 = -t50 * t20 + t61 * t22;
t9 = t61 * t20 + t50 * t22;
t4 = -t20 * pkin(10) + t6;
t3 = -t41 * pkin(5) - t22 * pkin(10) + t5;
t2 = t50 * t3 + t61 * t4;
t1 = t61 * t3 - t50 * t4;
t7 = [0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t49 * t64, -t48 * t64, 0, qJD(2) ^ 2 / 0.2e1 + (t48 ^ 2 / 0.2e1 + t49 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t55, t53 ^ 2 * t47, t53 * t55 * t54, t53 * t58, t54 ^ 2 * t47, t54 * t58, qJD(3) ^ 2 / 0.2e1, t29 * qJD(3) - t40 * t59, -t30 * qJD(3) + t40 * t60 (-t29 * t53 + t30 * t54) * qJD(1), t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1, t36 ^ 2 / 0.2e1, -t36 * t34, -t36 * t43, t34 ^ 2 / 0.2e1, t34 * t43, t43 ^ 2 / 0.2e1, -t16 * t43 + t26 * t34, t17 * t43 + t26 * t36, -t16 * t36 - t17 * t34, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, -t22 * t41, t20 ^ 2 / 0.2e1, t20 * t41, t41 ^ 2 / 0.2e1, t19 * t20 - t5 * t41, t19 * t22 + t6 * t41, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t9, -t11 * t37, t9 ^ 2 / 0.2e1, t9 * t37, t37 ^ 2 / 0.2e1, -t1 * t37 + t13 * t9, t13 * t11 + t2 * t37, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg  = t7;
