% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:06
% EndTime: 2019-03-09 15:30:06
% DurationCPUTime: 0.21s
% Computational Cost: add. (458->61), mult. (1032->126), div. (0->0), fcn. (643->6), ass. (0->50)
t44 = sin(qJ(3));
t45 = sin(qJ(2));
t46 = cos(qJ(3));
t47 = cos(qJ(2));
t26 = (t44 * t47 + t45 * t46) * qJD(1);
t48 = qJD(1) ^ 2;
t66 = t48 / 0.2e1;
t65 = -pkin(4) - pkin(9);
t64 = -pkin(8) - pkin(7);
t63 = cos(qJ(6));
t57 = qJD(1) * t47;
t58 = qJD(1) * t45;
t24 = t44 * t58 - t46 * t57;
t62 = t26 * t24;
t38 = qJD(2) + qJD(3);
t61 = t26 * t38;
t60 = t38 * t24;
t59 = t47 * t48;
t30 = qJD(2) * pkin(2) + t64 * t58;
t31 = t64 * t57;
t13 = t44 * t30 - t46 * t31;
t32 = -qJD(1) * pkin(1) - pkin(2) * t57;
t19 = t24 ^ 2 / 0.2e1;
t20 = t26 ^ 2 / 0.2e1;
t35 = t38 ^ 2 / 0.2e1;
t56 = qJD(1) * qJD(2);
t11 = t38 * qJ(4) + t13;
t55 = t45 * t56;
t54 = t47 * t56;
t12 = t46 * t30 + t44 * t31;
t9 = t24 * pkin(3) - t26 * qJ(4) + t32;
t53 = qJD(4) - t12;
t52 = qJD(5) - t9;
t8 = -t24 * qJ(5) - t11;
t50 = -t26 * qJ(5) + t53;
t43 = sin(qJ(6));
t41 = t47 ^ 2;
t40 = t45 ^ 2;
t23 = qJD(6) + t26;
t16 = t63 * t24 - t43 * t38;
t14 = t43 * t24 + t63 * t38;
t10 = -t38 * pkin(3) + t53;
t7 = t38 * pkin(5) - t8;
t6 = -t24 * pkin(4) + t52;
t5 = (-pkin(3) - pkin(4)) * t38 + t50;
t4 = (-pkin(3) + t65) * t38 + t50;
t3 = t26 * pkin(5) + t65 * t24 + t52;
t2 = t43 * t3 + t63 * t4;
t1 = t63 * t3 - t43 * t4;
t15 = [0, 0, 0, 0, 0, t66, 0, 0, 0, 0, t40 * t66, t45 * t59, t55, t41 * t66, t54, qJD(2) ^ 2 / 0.2e1, pkin(1) * t59 - pkin(7) * t55, -t48 * pkin(1) * t45 - pkin(7) * t54 (t40 + t41) * t48 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t41 / 0.2e1 + t40 / 0.2e1) * pkin(7) ^ 2) * t48, t20, -t62, t61, t19, -t60, t35, t12 * t38 + t32 * t24, -t13 * t38 + t32 * t26, -t12 * t26 - t13 * t24, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t20, t61, t62, t35, t60, t19, -t10 * t38 + t9 * t24, t10 * t26 - t11 * t24, t11 * t38 - t9 * t26, t11 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t19, -t62, -t60, t20, t61, t35, t6 * t26 - t8 * t38, t6 * t24 + t5 * t38, -t8 * t24 - t5 * t26, t5 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t23, t14 ^ 2 / 0.2e1, -t14 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t7 * t14, t7 * t16 - t2 * t23, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t15;
