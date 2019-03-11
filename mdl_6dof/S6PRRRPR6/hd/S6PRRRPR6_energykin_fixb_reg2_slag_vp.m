% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:36:45
% EndTime: 2019-03-08 23:36:45
% DurationCPUTime: 0.16s
% Computational Cost: add. (420->58), mult. (937->132), div. (0->0), fcn. (636->10), ass. (0->53)
t49 = qJD(2) ^ 2;
t67 = t49 / 0.2e1;
t66 = pkin(4) + pkin(5);
t65 = cos(qJ(4));
t64 = cos(qJ(6));
t44 = sin(qJ(4));
t45 = sin(qJ(3));
t59 = qJD(2) * t45;
t27 = -t65 * qJD(3) + t44 * t59;
t29 = t44 * qJD(3) + t65 * t59;
t63 = t29 * t27;
t47 = cos(qJ(3));
t58 = t47 * qJD(2);
t37 = -qJD(4) + t58;
t62 = t37 * t27;
t46 = sin(qJ(2));
t41 = sin(pkin(6));
t61 = qJD(1) * t41;
t30 = qJD(2) * pkin(8) + t46 * t61;
t42 = cos(pkin(6));
t60 = qJD(1) * t42;
t20 = t47 * t30 + t45 * t60;
t17 = qJD(3) * pkin(9) + t20;
t48 = cos(qJ(2));
t55 = t48 * t61;
t21 = -t55 + (-pkin(3) * t47 - pkin(9) * t45 - pkin(2)) * qJD(2);
t9 = t65 * t17 + t44 * t21;
t19 = -t45 * t30 + t47 * t60;
t57 = t27 ^ 2 / 0.2e1;
t56 = qJD(2) * qJD(3);
t7 = -t37 * qJ(5) + t9;
t54 = qJD(2) * t61;
t53 = qJD(3) * pkin(3) + t19;
t8 = -t44 * t17 + t65 * t21;
t52 = qJD(5) - t8;
t51 = t29 * qJ(5) + t53;
t50 = qJD(1) ^ 2;
t43 = sin(qJ(6));
t34 = qJD(6) + t37;
t32 = t37 ^ 2 / 0.2e1;
t31 = -qJD(2) * pkin(2) - t55;
t24 = t29 ^ 2 / 0.2e1;
t22 = t29 * t37;
t13 = t43 * t27 + t64 * t29;
t11 = -t64 * t27 + t43 * t29;
t10 = t27 * pkin(4) - t51;
t6 = t37 * pkin(4) + t52;
t5 = -t66 * t27 + t51;
t4 = t27 * pkin(10) + t7;
t3 = -t29 * pkin(10) + t66 * t37 + t52;
t2 = t43 * t3 + t64 * t4;
t1 = t64 * t3 - t43 * t4;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t50 / 0.2e1, 0, 0, 0, 0, 0, t67, t48 * t54, -t46 * t54, 0 (t42 ^ 2 / 0.2e1 + (t46 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1) * t41 ^ 2) * t50, t45 ^ 2 * t67, t45 * t49 * t47, t45 * t56, t47 ^ 2 * t67, t47 * t56, qJD(3) ^ 2 / 0.2e1, t19 * qJD(3) - t31 * t58, -t20 * qJD(3) + t31 * t59 (-t19 * t45 + t20 * t47) * qJD(2), t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t24, -t63, -t22, t57, t62, t32, -t27 * t53 - t8 * t37, -t29 * t53 + t9 * t37, -t9 * t27 - t8 * t29, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t53 ^ 2 / 0.2e1, t24, -t22, t63, t32, -t62, t57, t10 * t27 + t6 * t37, -t7 * t27 + t6 * t29, -t10 * t29 - t7 * t37, t7 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t34, t11 ^ 2 / 0.2e1, -t11 * t34, t34 ^ 2 / 0.2e1, t1 * t34 + t5 * t11, t5 * t13 - t2 * t34, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
