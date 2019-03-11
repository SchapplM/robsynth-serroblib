% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:25
% EndTime: 2019-03-09 02:57:25
% DurationCPUTime: 0.16s
% Computational Cost: add. (381->53), mult. (802->113), div. (0->0), fcn. (456->6), ass. (0->46)
t41 = qJD(1) ^ 2;
t33 = t41 / 0.2e1;
t56 = pkin(4) + pkin(8);
t55 = cos(qJ(6));
t37 = cos(pkin(9));
t39 = sin(qJ(3));
t40 = cos(qJ(3));
t52 = sin(pkin(9));
t21 = (-t37 * t39 - t52 * t40) * qJD(1);
t51 = qJD(1) * t39;
t23 = t37 * t40 * qJD(1) - t52 * t51;
t54 = t23 * t21;
t27 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t44 = -qJ(4) * qJD(1) + t27;
t17 = qJD(3) * pkin(3) + t44 * t40;
t19 = t44 * t39;
t9 = t52 * t17 + t37 * t19;
t53 = t41 * qJ(2);
t50 = qJD(3) * t21;
t49 = qJD(3) * t27;
t48 = t23 * qJD(3);
t47 = t21 ^ 2 / 0.2e1;
t46 = t23 ^ 2 / 0.2e1;
t45 = qJD(1) * qJD(3);
t25 = pkin(3) * t51 + qJD(1) * qJ(2) + qJD(4);
t8 = t37 * t17 - t52 * t19;
t43 = qJD(5) - t8;
t7 = -qJD(3) * qJ(5) - t9;
t42 = -t23 * qJ(5) + t25;
t38 = sin(qJ(6));
t36 = t40 ^ 2;
t35 = t39 ^ 2;
t32 = qJD(3) ^ 2 / 0.2e1;
t31 = qJ(2) ^ 2 * t33;
t30 = -pkin(1) * qJD(1) + qJD(2);
t20 = qJD(6) + t23;
t13 = t55 * qJD(3) - t38 * t21;
t11 = t38 * qJD(3) + t55 * t21;
t10 = -t21 * pkin(4) + t42;
t6 = -qJD(3) * pkin(4) + t43;
t5 = -t56 * t21 + t42;
t4 = t21 * pkin(5) - t7;
t3 = t23 * pkin(5) - t56 * qJD(3) + t43;
t2 = t38 * t3 + t55 * t5;
t1 = t55 * t3 - t38 * t5;
t12 = [0, 0, 0, 0, 0, t33, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, t30 * qJD(1), t53, t31 + t30 ^ 2 / 0.2e1, t36 * t33, -t40 * t41 * t39, t40 * t45, t35 * t33, -t39 * t45, t32, t39 * t53 + t40 * t49, -t39 * t49 + t40 * t53 (-t35 - t36) * t27 * qJD(1), t31 + (t35 / 0.2e1 + t36 / 0.2e1) * t27 ^ 2, t46, t54, t48, t47, t50, t32, t8 * qJD(3) - t25 * t21, -t9 * qJD(3) + t25 * t23, t9 * t21 - t8 * t23, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t32, -t48, -t50, t46, t54, t47, -t7 * t21 + t6 * t23, t6 * qJD(3) + t10 * t21, -t7 * qJD(3) - t10 * t23, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t20, t11 ^ 2 / 0.2e1, -t11 * t20, t20 ^ 2 / 0.2e1, t1 * t20 + t4 * t11, t4 * t13 - t2 * t20, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t12;
