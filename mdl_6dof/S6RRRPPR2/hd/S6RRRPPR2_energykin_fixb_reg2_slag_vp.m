% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:42
% EndTime: 2019-03-09 15:26:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (759->65), mult. (1851->141), div. (0->0), fcn. (1327->8), ass. (0->55)
t49 = qJD(1) ^ 2;
t68 = t49 / 0.2e1;
t67 = pkin(4) + pkin(9);
t66 = -pkin(8) - pkin(7);
t65 = cos(qJ(3));
t64 = cos(qJ(6));
t46 = sin(qJ(3));
t48 = cos(qJ(2));
t58 = qJD(1) * t48;
t47 = sin(qJ(2));
t59 = qJD(1) * t47;
t31 = t46 * t59 - t65 * t58;
t33 = (t46 * t48 + t65 * t47) * qJD(1);
t43 = sin(pkin(10));
t44 = cos(pkin(10));
t21 = t44 * t31 + t43 * t33;
t23 = -t43 * t31 + t44 * t33;
t63 = t23 * t21;
t40 = qJD(2) + qJD(3);
t62 = t23 * t40;
t61 = t40 * t21;
t60 = t48 * t49;
t35 = qJD(2) * pkin(2) + t66 * t59;
t36 = t66 * t58;
t25 = t65 * t35 + t46 * t36;
t12 = t40 * pkin(3) - t33 * qJ(4) + t25;
t26 = t46 * t35 - t65 * t36;
t15 = -t31 * qJ(4) + t26;
t9 = t43 * t12 + t44 * t15;
t57 = t21 ^ 2 / 0.2e1;
t56 = t23 ^ 2 / 0.2e1;
t55 = qJD(1) * qJD(2);
t54 = t47 * t55;
t53 = t48 * t55;
t8 = t44 * t12 - t43 * t15;
t7 = -t40 * qJ(5) - t9;
t52 = qJD(5) - t8;
t37 = (-pkin(2) * t48 - pkin(1)) * qJD(1);
t27 = t31 * pkin(3) + qJD(4) + t37;
t51 = -t23 * qJ(5) + t27;
t45 = sin(qJ(6));
t42 = t48 ^ 2;
t41 = t47 ^ 2;
t39 = t40 ^ 2 / 0.2e1;
t20 = qJD(6) + t23;
t18 = t45 * t21 + t64 * t40;
t16 = -t64 * t21 + t45 * t40;
t10 = t21 * pkin(4) + t51;
t6 = -t40 * pkin(4) + t52;
t5 = t67 * t21 + t51;
t4 = -t21 * pkin(5) - t7;
t3 = t23 * pkin(5) - t67 * t40 + t52;
t2 = t45 * t3 + t64 * t5;
t1 = t64 * t3 - t45 * t5;
t11 = [0, 0, 0, 0, 0, t68, 0, 0, 0, 0, t41 * t68, t47 * t60, t54, t42 * t68, t53, qJD(2) ^ 2 / 0.2e1, pkin(1) * t60 - pkin(7) * t54, -t49 * pkin(1) * t47 - pkin(7) * t53 (t41 + t42) * t49 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * pkin(7) ^ 2) * t49, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * t40, t31 ^ 2 / 0.2e1, -t31 * t40, t39, t25 * t40 + t37 * t31, -t26 * t40 + t37 * t33, -t25 * t33 - t26 * t31, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t56, -t63, t62, t57, -t61, t39, t27 * t21 + t8 * t40, t27 * t23 - t9 * t40, -t9 * t21 - t8 * t23, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t39, -t62, t61, t56, -t63, t57, t7 * t21 + t6 * t23, -t10 * t21 + t6 * t40, -t10 * t23 - t7 * t40, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t20, t16 ^ 2 / 0.2e1, -t16 * t20, t20 ^ 2 / 0.2e1, t1 * t20 + t4 * t16, t4 * t18 - t2 * t20, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
