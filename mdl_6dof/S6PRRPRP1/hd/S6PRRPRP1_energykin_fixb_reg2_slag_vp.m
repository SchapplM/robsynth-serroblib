% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:33
% EndTime: 2019-03-08 21:26:33
% DurationCPUTime: 0.15s
% Computational Cost: add. (462->58), mult. (1138->131), div. (0->0), fcn. (822->10), ass. (0->52)
t52 = qJD(2) ^ 2;
t64 = t52 / 0.2e1;
t50 = cos(qJ(3));
t51 = cos(qJ(2));
t45 = sin(pkin(6));
t61 = qJD(1) * t45;
t55 = t51 * t61;
t30 = -t55 + qJD(4) + (-pkin(3) * t50 - pkin(2)) * qJD(2);
t44 = sin(pkin(11));
t58 = qJD(2) * t50;
t48 = sin(qJ(3));
t59 = qJD(2) * t48;
t62 = cos(pkin(11));
t32 = t44 * t59 - t62 * t58;
t34 = (t44 * t50 + t62 * t48) * qJD(2);
t13 = t32 * pkin(4) - t34 * pkin(9) + t30;
t47 = sin(qJ(5));
t63 = cos(qJ(5));
t49 = sin(qJ(2));
t36 = qJD(2) * pkin(8) + t49 * t61;
t46 = cos(pkin(6));
t60 = qJD(1) * t46;
t40 = t50 * t60;
t57 = qJ(4) * qJD(2);
t20 = qJD(3) * pkin(3) + t40 + (-t36 - t57) * t48;
t28 = t50 * t36 + t48 * t60;
t21 = t50 * t57 + t28;
t10 = t44 * t20 + t62 * t21;
t8 = qJD(3) * pkin(9) + t10;
t4 = t47 * t13 + t63 * t8;
t56 = qJD(2) * qJD(3);
t3 = t63 * t13 - t47 * t8;
t54 = qJD(2) * t61;
t9 = t62 * t20 - t44 * t21;
t7 = -qJD(3) * pkin(4) - t9;
t53 = qJD(1) ^ 2;
t43 = qJD(3) ^ 2 / 0.2e1;
t37 = -qJD(2) * pkin(2) - t55;
t31 = qJD(5) + t32;
t29 = t31 ^ 2 / 0.2e1;
t27 = -t48 * t36 + t40;
t26 = t47 * qJD(3) + t63 * t34;
t24 = -t63 * qJD(3) + t47 * t34;
t23 = t26 ^ 2 / 0.2e1;
t22 = t24 ^ 2 / 0.2e1;
t16 = t26 * t31;
t15 = t24 * t31;
t14 = t26 * t24;
t5 = t24 * pkin(5) + qJD(6) + t7;
t2 = -t24 * qJ(6) + t4;
t1 = t31 * pkin(5) - t26 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t53 / 0.2e1, 0, 0, 0, 0, 0, t64, t51 * t54, -t49 * t54, 0 (t46 ^ 2 / 0.2e1 + (t49 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1) * t45 ^ 2) * t53, t48 ^ 2 * t64, t48 * t52 * t50, t48 * t56, t50 ^ 2 * t64, t50 * t56, t43, t27 * qJD(3) - t37 * t58, -t28 * qJD(3) + t37 * t59 (-t27 * t48 + t28 * t50) * qJD(2), t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * qJD(3), t32 ^ 2 / 0.2e1, -t32 * qJD(3), t43, t9 * qJD(3) + t30 * t32, -t10 * qJD(3) + t30 * t34, -t10 * t32 - t9 * t34, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t23, -t14, t16, t22, -t15, t29, t7 * t24 + t3 * t31, t7 * t26 - t4 * t31, -t4 * t24 - t3 * t26, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t23, -t14, t16, t22, -t15, t29, t1 * t31 + t5 * t24, -t2 * t31 + t5 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
