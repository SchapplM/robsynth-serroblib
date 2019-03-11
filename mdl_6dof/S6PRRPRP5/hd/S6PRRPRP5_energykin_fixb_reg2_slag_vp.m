% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:15
% EndTime: 2019-03-08 21:49:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (283->53), mult. (654->113), div. (0->0), fcn. (393->8), ass. (0->52)
t44 = qJD(2) ^ 2;
t61 = t44 / 0.2e1;
t60 = -pkin(3) - pkin(9);
t42 = cos(qJ(3));
t39 = sin(qJ(3));
t47 = -qJ(4) * t39 - pkin(2);
t43 = cos(qJ(2));
t36 = sin(pkin(6));
t57 = qJD(1) * t36;
t51 = t43 * t57;
t10 = -t51 + (t60 * t42 + t47) * qJD(2);
t38 = sin(qJ(5));
t41 = cos(qJ(5));
t40 = sin(qJ(2));
t23 = qJD(2) * pkin(8) + t40 * t57;
t37 = cos(pkin(6));
t56 = qJD(1) * t37;
t13 = -t39 * t23 + t42 * t56;
t46 = qJD(4) - t13;
t54 = t39 * qJD(2);
t7 = pkin(4) * t54 + t60 * qJD(3) + t46;
t4 = t41 * t10 + t38 * t7;
t55 = qJD(2) * t42;
t20 = t38 * qJD(3) + t41 * t55;
t22 = t41 * qJD(3) - t38 * t55;
t59 = t22 * t20;
t29 = qJD(5) + t54;
t58 = t29 * t20;
t14 = t42 * t23 + t39 * t56;
t53 = t20 ^ 2 / 0.2e1;
t52 = qJD(2) * qJD(3);
t12 = -qJD(3) * qJ(4) - t14;
t50 = qJD(2) * t57;
t49 = t39 * t52;
t48 = t42 * t52;
t8 = pkin(4) * t55 - t12;
t3 = -t38 * t10 + t41 * t7;
t45 = qJD(1) ^ 2;
t34 = qJD(3) ^ 2 / 0.2e1;
t31 = t42 ^ 2 * t61;
t30 = t39 ^ 2 * t61;
t28 = t39 * t44 * t42;
t25 = t29 ^ 2 / 0.2e1;
t24 = -qJD(2) * pkin(2) - t51;
t17 = t22 ^ 2 / 0.2e1;
t16 = t22 * t29;
t15 = -t51 + (-pkin(3) * t42 + t47) * qJD(2);
t11 = -qJD(3) * pkin(3) + t46;
t5 = t20 * pkin(5) - t22 * qJ(6) + t8;
t2 = t29 * qJ(6) + t4;
t1 = -t29 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t45 / 0.2e1, 0, 0, 0, 0, 0, t61, t43 * t50, -t40 * t50, 0 (t37 ^ 2 / 0.2e1 + (t40 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1) * t36 ^ 2) * t45, t30, t28, t49, t31, t48, t34, t13 * qJD(3) - t24 * t55, -t14 * qJD(3) + t24 * t54 (-t13 * t39 + t14 * t42) * qJD(2), t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t34, -t49, -t48, t30, t28, t31 (t11 * t39 - t12 * t42) * qJD(2), t11 * qJD(3) + t15 * t55, -t12 * qJD(3) - t15 * t54, t15 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t17, -t59, t16, t53, -t58, t25, t8 * t20 + t3 * t29, t8 * t22 - t4 * t29, -t4 * t20 - t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t17, t16, t59, t25, t58, t53, -t1 * t29 + t5 * t20, t1 * t22 - t2 * t20, t2 * t29 - t5 * t22, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
