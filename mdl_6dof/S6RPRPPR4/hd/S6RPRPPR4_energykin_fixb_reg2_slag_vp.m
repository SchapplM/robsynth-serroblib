% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:22
% EndTime: 2019-03-09 02:48:22
% DurationCPUTime: 0.22s
% Computational Cost: add. (624->62), mult. (1599->133), div. (0->0), fcn. (1141->8), ass. (0->51)
t52 = qJD(1) ^ 2;
t67 = t52 / 0.2e1;
t66 = -pkin(4) - pkin(5);
t65 = cos(qJ(3));
t64 = cos(qJ(6));
t48 = sin(pkin(9));
t49 = cos(pkin(9));
t51 = sin(qJ(3));
t38 = (t65 * t48 + t49 * t51) * qJD(1);
t47 = sin(pkin(10));
t60 = cos(pkin(10));
t25 = -t60 * qJD(3) + t47 * t38;
t27 = t47 * qJD(3) + t60 * t38;
t63 = t27 * t25;
t58 = qJD(1) * t49;
t59 = qJD(1) * t48;
t36 = t51 * t59 - t65 * t58;
t62 = t36 * t25;
t61 = pkin(7) + qJ(2);
t41 = qJD(2) + (-pkin(2) * t49 - pkin(1)) * qJD(1);
t15 = t36 * pkin(3) - t38 * qJ(4) + t41;
t39 = t61 * t59;
t40 = t61 * t58;
t22 = -t51 * t39 + t65 * t40;
t20 = qJD(3) * qJ(4) + t22;
t10 = t47 * t15 + t60 * t20;
t21 = -t65 * t39 - t51 * t40;
t57 = t25 ^ 2 / 0.2e1;
t28 = t36 ^ 2 / 0.2e1;
t7 = t36 * qJ(5) + t10;
t9 = t60 * t15 - t47 * t20;
t56 = qJD(3) * pkin(3) - qJD(4) + t21;
t55 = qJD(5) - t9;
t54 = t27 * qJ(5) + t56;
t50 = sin(qJ(6));
t46 = t49 ^ 2;
t45 = t48 ^ 2;
t43 = -qJD(1) * pkin(1) + qJD(2);
t30 = -qJD(6) + t36;
t24 = t27 ^ 2 / 0.2e1;
t16 = t27 * t36;
t13 = t50 * t25 + t64 * t27;
t11 = -t64 * t25 + t50 * t27;
t8 = t25 * pkin(4) - t54;
t6 = -t36 * pkin(4) + t55;
t5 = t66 * t25 + t54;
t4 = t25 * pkin(8) + t7;
t3 = -t27 * pkin(8) + t66 * t36 + t55;
t2 = t50 * t3 + t64 * t4;
t1 = t64 * t3 - t50 * t4;
t12 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t45 * t67, t48 * t52 * t49, 0, t46 * t67, 0, 0, -t43 * t58, t43 * t59 (t45 + t46) * t52 * qJ(2), t43 ^ 2 / 0.2e1 + (t46 / 0.2e1 + t45 / 0.2e1) * qJ(2) ^ 2 * t52, t38 ^ 2 / 0.2e1, -t38 * t36, t38 * qJD(3), t28, -t36 * qJD(3), qJD(3) ^ 2 / 0.2e1, t21 * qJD(3) + t41 * t36, -t22 * qJD(3) + t41 * t38, -t21 * t38 - t22 * t36, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1, t24, -t63, t16, t57, -t62, t28, -t25 * t56 + t9 * t36, -t10 * t36 - t27 * t56, -t10 * t25 - t9 * t27, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t56 ^ 2 / 0.2e1, t24, t16, t63, t28, t62, t57, t8 * t25 - t6 * t36, -t7 * t25 + t6 * t27, -t8 * t27 + t7 * t36, t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, -t13 * t30, t11 ^ 2 / 0.2e1, t11 * t30, t30 ^ 2 / 0.2e1, -t1 * t30 + t5 * t11, t5 * t13 + t2 * t30, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
