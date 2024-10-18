% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR15_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:17
% EndTime: 2024-09-27 22:26:17
% DurationCPUTime: 0.03s
% Computational Cost: add. (996->64), mult. (2978->152), div. (0->0), fcn. (2379->12), ass. (0->54)
t54 = sin(qJ(3));
t56 = cos(qJ(2));
t49 = sin(pkin(5));
t67 = qJD(1) * t49;
t62 = t56 * t67;
t55 = sin(qJ(2));
t63 = t55 * t67;
t71 = cos(qJ(3));
t32 = t54 * t63 - t71 * t62;
t34 = (t54 * t56 + t71 * t55) * t67;
t53 = sin(qJ(4));
t70 = cos(qJ(4));
t24 = -t53 * t32 + t70 * t34;
t72 = pkin(11) * t24;
t69 = cos(qJ(5));
t57 = qJD(1) ^ 2;
t68 = t49 ^ 2 * t57;
t66 = cos(pkin(5)) * qJD(1);
t64 = pkin(1) * t66;
t45 = t56 * t64;
t46 = qJD(2) + t66;
t28 = t46 * pkin(2) + t45 + (-pkin(8) - pkin(9)) * t63;
t36 = pkin(8) * t62 + t55 * t64;
t30 = pkin(9) * t62 + t36;
t19 = t71 * t28 - t54 * t30;
t42 = qJD(3) + t46;
t14 = t42 * pkin(3) - t34 * pkin(10) + t19;
t20 = t54 * t28 + t71 * t30;
t16 = -t32 * pkin(10) + t20;
t7 = t53 * t14 + t70 * t16;
t65 = t56 * t68;
t61 = t68 / 0.2e1;
t48 = sin(pkin(6));
t60 = t48 * t69;
t50 = cos(pkin(6));
t59 = t50 * t69;
t6 = t70 * t14 - t53 * t16;
t22 = t70 * t32 + t53 * t34;
t41 = qJD(4) + t42;
t58 = -t22 * t50 + t41 * t48;
t38 = (-pkin(2) * t56 - pkin(1)) * t67;
t27 = t32 * pkin(3) + t38;
t52 = sin(qJ(5));
t35 = -pkin(8) * t63 + t45;
t17 = -t48 * t22 - t50 * t41 - qJD(5);
t11 = t69 * t24 + t58 * t52;
t9 = t22 * t59 + t52 * t24 - t41 * t60;
t8 = t22 * pkin(4) - t48 * t72 + t27;
t5 = t41 * pkin(4) - t50 * t72 + t6;
t4 = t58 * pkin(11) + t7;
t3 = -t48 * t5 + t50 * t8;
t2 = t69 * t4 + (t48 * t8 + t5 * t50) * t52;
t1 = -t52 * t4 + t5 * t59 + t8 * t60;
t10 = [0, 0, 0, 0, 0, t57 / 0.2e1, 0, 0, 0, 0, t55 ^ 2 * t61, t55 * t65, t46 * t63, t56 ^ 2 * t61, t46 * t62, t46 ^ 2 / 0.2e1, pkin(1) * t65 + t35 * t46, -pkin(1) * t55 * t68 - t36 * t46, (-t35 * t55 + t36 * t56) * t67, t36 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t61, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t42, t32 ^ 2 / 0.2e1, -t32 * t42, t42 ^ 2 / 0.2e1, t19 * t42 + t38 * t32, -t20 * t42 + t38 * t34, -t19 * t34 - t20 * t32, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, t24 * t41, t22 ^ 2 / 0.2e1, -t22 * t41, t41 ^ 2 / 0.2e1, t27 * t22 + t6 * t41, t27 * t24 - t7 * t41, -t7 * t22 - t6 * t24, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t9, -t11 * t17, t9 ^ 2 / 0.2e1, t9 * t17, t17 ^ 2 / 0.2e1, -t1 * t17 + t3 * t9, t3 * t11 + t2 * t17, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t10;
