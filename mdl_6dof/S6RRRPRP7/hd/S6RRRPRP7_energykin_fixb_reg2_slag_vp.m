% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:11:06
% EndTime: 2019-03-09 17:11:06
% DurationCPUTime: 0.25s
% Computational Cost: add. (1255->66), mult. (2994->147), div. (0->0), fcn. (2344->10), ass. (0->55)
t59 = cos(qJ(2));
t58 = sin(qJ(2));
t54 = sin(pkin(6));
t68 = qJD(1) * t54;
t63 = t58 * t68;
t67 = cos(pkin(6)) * qJD(1);
t64 = pkin(1) * t67;
t41 = -pkin(8) * t63 + t59 * t64;
t51 = qJD(2) + t67;
t34 = -t51 * pkin(2) - t41;
t57 = sin(qJ(3));
t74 = cos(qJ(3));
t38 = -t74 * t51 + t57 * t63;
t27 = t38 * pkin(3) + qJD(4) + t34;
t40 = t57 * t51 + t74 * t63;
t53 = sin(pkin(11));
t69 = cos(pkin(11));
t28 = t69 * t38 + t53 * t40;
t30 = -t53 * t38 + t69 * t40;
t12 = t28 * pkin(4) - t30 * pkin(10) + t27;
t56 = sin(qJ(5));
t73 = cos(qJ(5));
t62 = t59 * t68;
t42 = pkin(8) * t62 + t58 * t64;
t35 = t51 * pkin(9) + t42;
t37 = (-pkin(2) * t59 - pkin(9) * t58 - pkin(1)) * t68;
t23 = -t57 * t35 + t74 * t37;
t46 = -qJD(3) + t62;
t15 = -t46 * pkin(3) - t40 * qJ(4) + t23;
t24 = t74 * t35 + t57 * t37;
t18 = -t38 * qJ(4) + t24;
t10 = t53 * t15 + t69 * t18;
t8 = -t46 * pkin(10) + t10;
t4 = t56 * t12 + t73 * t8;
t20 = t56 * t30 + t73 * t46;
t22 = t73 * t30 - t56 * t46;
t72 = t22 * t20;
t26 = qJD(5) + t28;
t71 = t26 * t20;
t60 = qJD(1) ^ 2;
t70 = t54 ^ 2 * t60;
t66 = t20 ^ 2 / 0.2e1;
t65 = t59 * t70;
t61 = t70 / 0.2e1;
t9 = t69 * t15 - t53 * t18;
t3 = t73 * t12 - t56 * t8;
t7 = t46 * pkin(4) - t9;
t44 = t46 ^ 2 / 0.2e1;
t25 = t26 ^ 2 / 0.2e1;
t19 = t22 ^ 2 / 0.2e1;
t13 = t22 * t26;
t5 = t20 * pkin(5) - t22 * qJ(6) + t7;
t2 = t26 * qJ(6) + t4;
t1 = -t26 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t60 / 0.2e1, 0, 0, 0, 0, t58 ^ 2 * t61, t58 * t65, t51 * t63, t59 ^ 2 * t61, t51 * t62, t51 ^ 2 / 0.2e1, pkin(1) * t65 + t41 * t51, -pkin(1) * t58 * t70 - t42 * t51 (-t41 * t58 + t42 * t59) * t68, t42 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t61, t40 ^ 2 / 0.2e1, -t40 * t38, -t40 * t46, t38 ^ 2 / 0.2e1, t38 * t46, t44, -t23 * t46 + t34 * t38, t24 * t46 + t34 * t40, -t23 * t40 - t24 * t38, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, -t30 * t46, t28 ^ 2 / 0.2e1, t28 * t46, t44, t27 * t28 - t9 * t46, t10 * t46 + t27 * t30, -t10 * t28 - t9 * t30, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t19, -t72, t13, t66, -t71, t25, t7 * t20 + t3 * t26, t7 * t22 - t4 * t26, -t4 * t20 - t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t19, t13, t72, t25, t71, t66, -t1 * t26 + t5 * t20, t1 * t22 - t2 * t20, t2 * t26 - t5 * t22, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
