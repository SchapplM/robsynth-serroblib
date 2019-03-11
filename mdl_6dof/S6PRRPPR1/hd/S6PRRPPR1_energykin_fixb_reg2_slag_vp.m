% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:11
% EndTime: 2019-03-08 21:02:11
% DurationCPUTime: 0.17s
% Computational Cost: add. (619->62), mult. (1509->148), div. (0->0), fcn. (1124->12), ass. (0->53)
t53 = qJD(2) ^ 2;
t67 = t53 / 0.2e1;
t50 = sin(qJ(2));
t46 = sin(pkin(6));
t63 = qJD(1) * t46;
t36 = qJD(2) * pkin(8) + t50 * t63;
t51 = cos(qJ(3));
t47 = cos(pkin(6));
t62 = qJD(1) * t47;
t40 = t51 * t62;
t49 = sin(qJ(3));
t58 = qJ(4) * qJD(2);
t22 = qJD(3) * pkin(3) + t40 + (-t36 - t58) * t49;
t29 = t51 * t36 + t49 * t62;
t23 = t51 * t58 + t29;
t45 = sin(pkin(11));
t65 = cos(pkin(11));
t12 = t45 * t22 + t65 * t23;
t10 = qJD(3) * qJ(5) + t12;
t52 = cos(qJ(2));
t56 = t52 * t63;
t30 = -t56 + qJD(4) + (-pkin(3) * t51 - pkin(2)) * qJD(2);
t60 = qJD(2) * t51;
t61 = qJD(2) * t49;
t32 = t45 * t61 - t65 * t60;
t34 = (t45 * t51 + t65 * t49) * qJD(2);
t18 = t32 * pkin(4) - t34 * qJ(5) + t30;
t44 = sin(pkin(12));
t64 = cos(pkin(12));
t6 = t64 * t10 + t44 * t18;
t66 = cos(qJ(6));
t59 = t32 ^ 2 / 0.2e1;
t57 = qJD(2) * qJD(3);
t55 = qJD(2) * t63;
t5 = -t44 * t10 + t64 * t18;
t11 = t65 * t22 - t45 * t23;
t9 = -qJD(3) * pkin(4) + qJD(5) - t11;
t54 = qJD(1) ^ 2;
t48 = sin(qJ(6));
t43 = qJD(3) ^ 2 / 0.2e1;
t37 = -qJD(2) * pkin(2) - t56;
t31 = qJD(6) + t32;
t28 = -t49 * t36 + t40;
t27 = t44 * qJD(3) + t64 * t34;
t25 = -t64 * qJD(3) + t44 * t34;
t15 = -t48 * t25 + t66 * t27;
t13 = t66 * t25 + t48 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(9) + t6;
t3 = t32 * pkin(5) - t27 * pkin(9) + t5;
t2 = t48 * t3 + t66 * t4;
t1 = t66 * t3 - t48 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t54 / 0.2e1, 0, 0, 0, 0, 0, t67, t52 * t55, -t50 * t55, 0 (t47 ^ 2 / 0.2e1 + (t50 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1) * t46 ^ 2) * t54, t49 ^ 2 * t67, t49 * t53 * t51, t49 * t57, t51 ^ 2 * t67, t51 * t57, t43, t28 * qJD(3) - t37 * t60, -t29 * qJD(3) + t37 * t61 (-t28 * t49 + t29 * t51) * qJD(2), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * qJD(3), t59, -t32 * qJD(3), t43, t11 * qJD(3) + t30 * t32, -t12 * qJD(3) + t30 * t34, -t11 * t34 - t12 * t32, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t32, t25 ^ 2 / 0.2e1, -t25 * t32, t59, t9 * t25 + t5 * t32, t9 * t27 - t6 * t32, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t31, t13 ^ 2 / 0.2e1, -t13 * t31, t31 ^ 2 / 0.2e1, t1 * t31 + t7 * t13, t7 * t15 - t2 * t31, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
