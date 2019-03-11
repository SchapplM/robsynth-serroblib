% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:38:28
% EndTime: 2019-03-09 21:38:28
% DurationCPUTime: 0.20s
% Computational Cost: add. (841->63), mult. (1980->133), div. (0->0), fcn. (1489->8), ass. (0->51)
t69 = -pkin(4) - pkin(5);
t68 = cos(qJ(3));
t67 = cos(qJ(4));
t63 = cos(pkin(6)) * qJD(1);
t47 = qJD(2) + t63;
t52 = sin(qJ(3));
t53 = sin(qJ(2));
t49 = sin(pkin(6));
t64 = qJD(1) * t49;
t60 = t53 * t64;
t37 = t52 * t47 + t68 * t60;
t54 = cos(qJ(2));
t59 = t54 * t64;
t42 = -qJD(3) + t59;
t51 = sin(qJ(4));
t21 = t51 * t37 + t67 * t42;
t23 = t67 * t37 - t51 * t42;
t9 = t23 * t21;
t35 = -t68 * t47 + t52 * t60;
t33 = qJD(4) + t35;
t16 = t23 * t33;
t66 = t33 * t21;
t55 = qJD(1) ^ 2;
t65 = t49 ^ 2 * t55;
t61 = pkin(1) * t63;
t38 = -pkin(8) * t60 + t54 * t61;
t28 = -t47 * pkin(2) - t38;
t11 = t35 * pkin(3) - t37 * pkin(10) + t28;
t39 = pkin(8) * t59 + t53 * t61;
t29 = t47 * pkin(9) + t39;
t32 = (-pkin(2) * t54 - pkin(9) * t53 - pkin(1)) * t64;
t18 = t68 * t29 + t52 * t32;
t15 = -t42 * pkin(10) + t18;
t8 = t51 * t11 + t67 * t15;
t17 = -t52 * t29 + t68 * t32;
t19 = t21 ^ 2 / 0.2e1;
t30 = t33 ^ 2 / 0.2e1;
t62 = t54 * t65;
t5 = t33 * qJ(5) + t8;
t58 = t65 / 0.2e1;
t14 = t42 * pkin(3) - t17;
t7 = t67 * t11 - t51 * t15;
t57 = qJD(5) - t7;
t56 = t23 * qJ(5) - t14;
t20 = t23 ^ 2 / 0.2e1;
t6 = t21 * pkin(4) - t56;
t4 = -t33 * pkin(4) + t57;
t3 = t69 * t21 + qJD(6) + t56;
t2 = t21 * qJ(6) + t5;
t1 = -t23 * qJ(6) + t69 * t33 + t57;
t10 = [0, 0, 0, 0, 0, t55 / 0.2e1, 0, 0, 0, 0, t53 ^ 2 * t58, t53 * t62, t47 * t60, t54 ^ 2 * t58, t47 * t59, t47 ^ 2 / 0.2e1, pkin(1) * t62 + t38 * t47, -pkin(1) * t53 * t65 - t39 * t47 (-t38 * t53 + t39 * t54) * t64, t39 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t58, t37 ^ 2 / 0.2e1, -t37 * t35, -t37 * t42, t35 ^ 2 / 0.2e1, t35 * t42, t42 ^ 2 / 0.2e1, -t17 * t42 + t28 * t35, t18 * t42 + t28 * t37, -t17 * t37 - t18 * t35, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t20, -t9, t16, t19, -t66, t30, t14 * t21 + t7 * t33, t14 * t23 - t8 * t33, -t8 * t21 - t7 * t23, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t20, t16, t9, t30, t66, t19, t6 * t21 - t4 * t33, -t5 * t21 + t4 * t23, -t6 * t23 + t5 * t33, t5 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t20, t9, -t16, t19, -t66, t30, -t1 * t33 - t3 * t21, t2 * t33 + t3 * t23, -t1 * t23 + t2 * t21, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t10;
