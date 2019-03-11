% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:03:58
% EndTime: 2019-03-09 22:03:59
% DurationCPUTime: 0.20s
% Computational Cost: add. (788->65), mult. (1851->143), div. (0->0), fcn. (1327->8), ass. (0->56)
t50 = qJD(1) ^ 2;
t69 = t50 / 0.2e1;
t68 = pkin(4) + pkin(10);
t67 = -pkin(8) - pkin(7);
t66 = cos(qJ(3));
t65 = cos(qJ(6));
t46 = sin(qJ(3));
t49 = cos(qJ(2));
t59 = qJD(1) * t49;
t47 = sin(qJ(2));
t60 = qJD(1) * t47;
t31 = t46 * t60 - t66 * t59;
t33 = (t46 * t49 + t66 * t47) * qJD(1);
t45 = sin(qJ(4));
t48 = cos(qJ(4));
t21 = t48 * t31 + t45 * t33;
t23 = -t45 * t31 + t48 * t33;
t64 = t23 * t21;
t41 = qJD(2) + qJD(3);
t40 = qJD(4) + t41;
t63 = t23 * t40;
t62 = t40 * t21;
t61 = t49 * t50;
t35 = qJD(2) * pkin(2) + t67 * t60;
t36 = t67 * t59;
t25 = t66 * t35 + t46 * t36;
t12 = t41 * pkin(3) - t33 * pkin(9) + t25;
t26 = t46 * t35 - t66 * t36;
t15 = -t31 * pkin(9) + t26;
t9 = t45 * t12 + t48 * t15;
t58 = t21 ^ 2 / 0.2e1;
t57 = t23 ^ 2 / 0.2e1;
t56 = qJD(1) * qJD(2);
t55 = t47 * t56;
t54 = t49 * t56;
t8 = t48 * t12 - t45 * t15;
t7 = -t40 * qJ(5) - t9;
t53 = qJD(5) - t8;
t37 = (-pkin(2) * t49 - pkin(1)) * qJD(1);
t27 = t31 * pkin(3) + t37;
t52 = -t23 * qJ(5) + t27;
t44 = sin(qJ(6));
t43 = t49 ^ 2;
t42 = t47 ^ 2;
t39 = t40 ^ 2 / 0.2e1;
t20 = qJD(6) + t23;
t18 = t44 * t21 + t65 * t40;
t16 = -t65 * t21 + t44 * t40;
t10 = t21 * pkin(4) + t52;
t6 = -t40 * pkin(4) + t53;
t5 = t68 * t21 + t52;
t4 = -t21 * pkin(5) - t7;
t3 = t23 * pkin(5) - t68 * t40 + t53;
t2 = t44 * t3 + t65 * t5;
t1 = t65 * t3 - t44 * t5;
t11 = [0, 0, 0, 0, 0, t69, 0, 0, 0, 0, t42 * t69, t47 * t61, t55, t43 * t69, t54, qJD(2) ^ 2 / 0.2e1, pkin(1) * t61 - pkin(7) * t55, -t50 * pkin(1) * t47 - pkin(7) * t54 (t42 + t43) * t50 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t43 / 0.2e1 + t42 / 0.2e1) * pkin(7) ^ 2) * t50, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * t41, t31 ^ 2 / 0.2e1, -t31 * t41, t41 ^ 2 / 0.2e1, t25 * t41 + t37 * t31, -t26 * t41 + t37 * t33, -t25 * t33 - t26 * t31, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t57, -t64, t63, t58, -t62, t39, t27 * t21 + t8 * t40, t27 * t23 - t9 * t40, -t9 * t21 - t8 * t23, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t39, -t63, t62, t57, -t64, t58, t7 * t21 + t6 * t23, -t10 * t21 + t6 * t40, -t10 * t23 - t7 * t40, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t20, t16 ^ 2 / 0.2e1, -t16 * t20, t20 ^ 2 / 0.2e1, t1 * t20 + t4 * t16, t4 * t18 - t2 * t20, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
