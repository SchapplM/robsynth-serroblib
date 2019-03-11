% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:05:28
% EndTime: 2019-03-09 14:05:28
% DurationCPUTime: 0.24s
% Computational Cost: add. (1337->69), mult. (3065->162), div. (0->0), fcn. (2249->10), ass. (0->53)
t62 = qJD(1) ^ 2;
t74 = t62 / 0.2e1;
t73 = cos(qJ(4));
t72 = cos(qJ(5));
t71 = cos(qJ(6));
t61 = cos(qJ(2));
t70 = t61 * t62;
t60 = sin(qJ(2));
t38 = (-pkin(2) * t61 - qJ(3) * t60 - pkin(1)) * qJD(1);
t67 = t61 * qJD(1);
t46 = pkin(7) * t67 + qJD(2) * qJ(3);
t56 = sin(pkin(11));
t69 = cos(pkin(11));
t32 = t69 * t38 - t56 * t46;
t68 = qJD(1) * t60;
t41 = t56 * qJD(2) + t69 * t68;
t25 = -pkin(3) * t67 - t41 * pkin(8) + t32;
t33 = t56 * t38 + t69 * t46;
t39 = -t69 * qJD(2) + t56 * t68;
t27 = -t39 * pkin(8) + t33;
t59 = sin(qJ(4));
t16 = t73 * t25 - t59 * t27;
t31 = -t59 * t39 + t73 * t41;
t49 = -qJD(4) + t67;
t12 = -t49 * pkin(4) - t31 * pkin(9) + t16;
t17 = t59 * t25 + t73 * t27;
t29 = t73 * t39 + t59 * t41;
t14 = -t29 * pkin(9) + t17;
t58 = sin(qJ(5));
t6 = t58 * t12 + t72 * t14;
t66 = qJD(1) * qJD(2);
t65 = t60 * t66;
t64 = t61 * t66;
t5 = t72 * t12 - t58 * t14;
t43 = -qJD(2) * pkin(2) + pkin(7) * t68 + qJD(3);
t47 = -qJD(5) + t49;
t34 = t39 * pkin(3) + t43;
t22 = t29 * pkin(4) + t34;
t57 = sin(qJ(6));
t55 = t61 ^ 2;
t54 = t60 ^ 2;
t51 = t55 * t74;
t44 = -qJD(6) + t47;
t21 = -t58 * t29 + t72 * t31;
t19 = t72 * t29 + t58 * t31;
t15 = t19 * pkin(5) + t22;
t9 = -t57 * t19 + t71 * t21;
t7 = t71 * t19 + t57 * t21;
t4 = -t19 * pkin(10) + t6;
t3 = -t47 * pkin(5) - t21 * pkin(10) + t5;
t2 = t57 * t3 + t71 * t4;
t1 = t71 * t3 - t57 * t4;
t8 = [0, 0, 0, 0, 0, t74, 0, 0, 0, 0, t54 * t74, t60 * t70, t65, t51, t64, qJD(2) ^ 2 / 0.2e1, pkin(1) * t70 - pkin(7) * t65, -t62 * pkin(1) * t60 - pkin(7) * t64 (t54 + t55) * t62 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t55 / 0.2e1 + t54 / 0.2e1) * pkin(7) ^ 2) * t62, t41 ^ 2 / 0.2e1, -t41 * t39, -t41 * t67, t39 ^ 2 / 0.2e1, t39 * t67, t51, -t32 * t67 + t43 * t39, t33 * t67 + t43 * t41, -t32 * t41 - t33 * t39, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t49, t29 ^ 2 / 0.2e1, t29 * t49, t49 ^ 2 / 0.2e1, -t16 * t49 + t34 * t29, t17 * t49 + t34 * t31, -t16 * t31 - t17 * t29, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, -t21 * t47, t19 ^ 2 / 0.2e1, t19 * t47, t47 ^ 2 / 0.2e1, t22 * t19 - t5 * t47, t22 * t21 + t6 * t47, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, -t9 * t44, t7 ^ 2 / 0.2e1, t7 * t44, t44 ^ 2 / 0.2e1, -t1 * t44 + t15 * t7, t15 * t9 + t2 * t44, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1;];
T_reg  = t8;
