% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPR4
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:35:59
% EndTime: 2019-03-09 15:35:59
% DurationCPUTime: 0.21s
% Computational Cost: add. (784->66), mult. (1732->142), div. (0->0), fcn. (1181->8), ass. (0->54)
t52 = qJD(1) ^ 2;
t70 = t52 / 0.2e1;
t69 = pkin(4) + pkin(5);
t68 = cos(qJ(3));
t67 = cos(qJ(6));
t49 = sin(qJ(3));
t50 = sin(qJ(2));
t62 = qJD(1) * t50;
t32 = -t68 * qJD(2) + t49 * t62;
t34 = t49 * qJD(2) + t68 * t62;
t47 = sin(pkin(10));
t63 = cos(pkin(10));
t22 = t63 * t32 + t47 * t34;
t24 = -t47 * t32 + t63 * t34;
t66 = t24 * t22;
t51 = cos(qJ(2));
t61 = t51 * qJD(1);
t41 = -qJD(3) + t61;
t65 = t41 * t22;
t64 = t51 * t52;
t31 = (-pkin(2) * t51 - pkin(8) * t50 - pkin(1)) * qJD(1);
t37 = pkin(7) * t61 + qJD(2) * pkin(8);
t25 = t68 * t31 - t49 * t37;
t15 = -t41 * pkin(3) - t34 * qJ(4) + t25;
t26 = t49 * t31 + t68 * t37;
t18 = -t32 * qJ(4) + t26;
t9 = t47 * t15 + t63 * t18;
t60 = t22 ^ 2 / 0.2e1;
t59 = qJD(1) * qJD(2);
t7 = -t41 * qJ(5) + t9;
t58 = t50 * t59;
t57 = t51 * t59;
t36 = -qJD(2) * pkin(2) + pkin(7) * t62;
t8 = t63 * t15 - t47 * t18;
t56 = qJD(5) - t8;
t55 = -t32 * pkin(3) - qJD(4) - t36;
t54 = t24 * qJ(5) + t55;
t48 = sin(qJ(6));
t46 = t51 ^ 2;
t45 = t50 ^ 2;
t40 = qJD(6) + t41;
t38 = t41 ^ 2 / 0.2e1;
t21 = t24 ^ 2 / 0.2e1;
t19 = t24 * t41;
t13 = t48 * t22 + t67 * t24;
t11 = -t67 * t22 + t48 * t24;
t10 = t22 * pkin(4) - t54;
t6 = t41 * pkin(4) + t56;
t5 = -t69 * t22 + t54;
t4 = t22 * pkin(9) + t7;
t3 = -t24 * pkin(9) + t69 * t41 + t56;
t2 = t48 * t3 + t67 * t4;
t1 = t67 * t3 - t48 * t4;
t12 = [0, 0, 0, 0, 0, t70, 0, 0, 0, 0, t45 * t70, t50 * t64, t58, t46 * t70, t57, qJD(2) ^ 2 / 0.2e1, pkin(1) * t64 - pkin(7) * t58, -t52 * pkin(1) * t50 - pkin(7) * t57 (t45 + t46) * t52 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t46 / 0.2e1 + t45 / 0.2e1) * pkin(7) ^ 2) * t52, t34 ^ 2 / 0.2e1, -t34 * t32, -t34 * t41, t32 ^ 2 / 0.2e1, t32 * t41, t38, -t25 * t41 + t36 * t32, t26 * t41 + t36 * t34, -t25 * t34 - t26 * t32, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t21, -t66, -t19, t60, t65, t38, -t22 * t55 - t8 * t41, -t24 * t55 + t9 * t41, -t9 * t22 - t8 * t24, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t55 ^ 2 / 0.2e1, t21, -t19, t66, t38, -t65, t60, t10 * t22 + t6 * t41, -t7 * t22 + t6 * t24, -t10 * t24 - t7 * t41, t7 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t40, t11 ^ 2 / 0.2e1, -t11 * t40, t40 ^ 2 / 0.2e1, t1 * t40 + t5 * t11, t5 * t13 - t2 * t40, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
