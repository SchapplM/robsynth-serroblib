% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:34
% EndTime: 2019-03-09 03:49:34
% DurationCPUTime: 0.18s
% Computational Cost: add. (590->61), mult. (1490->134), div. (0->0), fcn. (1041->8), ass. (0->50)
t55 = qJD(1) ^ 2;
t67 = t55 / 0.2e1;
t49 = sin(pkin(10));
t50 = cos(pkin(10));
t54 = sin(qJ(3));
t66 = cos(qJ(3));
t35 = (t66 * t49 + t50 * t54) * qJD(1);
t61 = qJD(1) * t49;
t62 = pkin(7) + qJ(2);
t36 = t62 * t61;
t60 = qJD(1) * t50;
t37 = t62 * t60;
t23 = -t66 * t36 - t54 * t37;
t57 = qJD(4) - t23;
t10 = -t35 * pkin(8) + (-pkin(3) - pkin(4)) * qJD(3) + t57;
t24 = -t54 * t36 + t66 * t37;
t22 = qJD(3) * qJ(4) + t24;
t33 = t54 * t61 - t66 * t60;
t12 = t33 * pkin(8) + t22;
t53 = sin(qJ(5));
t65 = cos(qJ(5));
t7 = t53 * t10 + t65 * t12;
t64 = cos(qJ(6));
t63 = t35 * t33;
t59 = qJD(3) * t33;
t58 = t33 ^ 2 / 0.2e1;
t42 = -qJD(1) * pkin(1) + qJD(2);
t38 = -pkin(2) * t60 + t42;
t18 = -t65 * t33 + t53 * t35;
t16 = t33 * pkin(3) - t35 * qJ(4) + t38;
t6 = t65 * t10 - t53 * t12;
t8 = -t33 * pkin(4) - t16;
t52 = sin(qJ(6));
t47 = qJD(3) ^ 2 / 0.2e1;
t45 = qJD(3) - qJD(5);
t44 = t50 ^ 2;
t43 = t49 ^ 2;
t27 = t35 * qJD(3);
t26 = t35 ^ 2 / 0.2e1;
t21 = -qJD(3) * pkin(3) + t57;
t20 = t53 * t33 + t65 * t35;
t17 = qJD(6) + t18;
t15 = t64 * t20 - t52 * t45;
t13 = t52 * t20 + t64 * t45;
t5 = -t45 * pkin(9) + t7;
t4 = t45 * pkin(5) - t6;
t3 = t18 * pkin(5) - t20 * pkin(9) + t8;
t2 = t52 * t3 + t64 * t5;
t1 = t64 * t3 - t52 * t5;
t9 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t43 * t67, t49 * t55 * t50, 0, t44 * t67, 0, 0, -t42 * t60, t42 * t61 (t43 + t44) * t55 * qJ(2), t42 ^ 2 / 0.2e1 + (t44 / 0.2e1 + t43 / 0.2e1) * qJ(2) ^ 2 * t55, t26, -t63, t27, t58, -t59, t47, t23 * qJD(3) + t38 * t33, -t24 * qJD(3) + t38 * t35, -t23 * t35 - t24 * t33, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t26, t27, t63, t47, t59, t58, -t21 * qJD(3) + t16 * t33, t21 * t35 - t22 * t33, t22 * qJD(3) - t16 * t35, t22 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, -t20 * t45, t18 ^ 2 / 0.2e1, t18 * t45, t45 ^ 2 / 0.2e1, t8 * t18 - t6 * t45, t8 * t20 + t7 * t45, -t7 * t18 - t6 * t20, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t17, t13 ^ 2 / 0.2e1, -t13 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t4 * t13, t4 * t15 - t2 * t17, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t9;
