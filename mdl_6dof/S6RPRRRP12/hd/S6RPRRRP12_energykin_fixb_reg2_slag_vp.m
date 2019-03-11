% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:50:37
% EndTime: 2019-03-09 06:50:37
% DurationCPUTime: 0.32s
% Computational Cost: add. (1671->76), mult. (5402->177), div. (0->0), fcn. (4449->12), ass. (0->65)
t58 = sin(pkin(12));
t61 = cos(pkin(12));
t60 = sin(pkin(6));
t77 = qJD(1) * t60;
t70 = qJ(2) * t77;
t63 = cos(pkin(6));
t76 = qJD(1) * t63;
t73 = pkin(1) * t76;
t49 = t58 * t73 + t61 * t70;
t59 = sin(pkin(7));
t62 = cos(pkin(7));
t79 = t60 * t61;
t36 = (t59 * t63 + t62 * t79) * qJD(1) * pkin(9) + t49;
t54 = t61 * t73;
t81 = t58 * t60;
t38 = t54 + (pkin(2) * t63 + (-pkin(9) * t62 - qJ(2)) * t81) * qJD(1);
t44 = qJD(2) + (-pkin(9) * t58 * t59 - pkin(2) * t61 - pkin(1)) * t77;
t66 = sin(qJ(3));
t67 = cos(qJ(3));
t20 = -t66 * t36 + (t38 * t62 + t44 * t59) * t67;
t68 = qJD(1) ^ 2;
t87 = t68 / 0.2e1;
t72 = t61 * t77;
t46 = t59 * t72 - t62 * t76 - qJD(3);
t18 = t46 * pkin(3) - t20;
t78 = t62 * t66;
t80 = t59 * t66;
t41 = (t63 * t80 + (t58 * t67 + t61 * t78) * t60) * qJD(1);
t65 = sin(qJ(4));
t86 = cos(qJ(4));
t29 = t65 * t41 + t86 * t46;
t31 = t86 * t41 - t65 * t46;
t12 = t29 * pkin(4) - t31 * pkin(11) + t18;
t64 = sin(qJ(5));
t26 = -t59 * t38 + t62 * t44;
t39 = t66 * t58 * t77 + (-t59 * t76 - t62 * t72) * t67;
t15 = t39 * pkin(3) - t41 * pkin(10) + t26;
t21 = t67 * t36 + t38 * t78 + t44 * t80;
t19 = -t46 * pkin(10) + t21;
t10 = t65 * t15 + t86 * t19;
t37 = qJD(4) + t39;
t8 = t37 * pkin(11) + t10;
t85 = cos(qJ(5));
t4 = t64 * t12 + t85 * t8;
t23 = t64 * t31 - t85 * t37;
t25 = t85 * t31 + t64 * t37;
t84 = t25 * t23;
t28 = qJD(5) + t29;
t83 = t28 * t23;
t82 = t60 ^ 2 * t68;
t75 = t23 ^ 2 / 0.2e1;
t74 = t60 * t63 * t68;
t71 = t82 / 0.2e1;
t9 = t86 * t15 - t65 * t19;
t3 = t85 * t12 - t64 * t8;
t7 = -t37 * pkin(4) - t9;
t55 = -pkin(1) * t77 + qJD(2);
t48 = -t58 * t70 + t54;
t27 = t28 ^ 2 / 0.2e1;
t22 = t25 ^ 2 / 0.2e1;
t13 = t25 * t28;
t5 = t23 * pkin(5) - t25 * qJ(6) + t7;
t2 = t28 * qJ(6) + t4;
t1 = -t28 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t87, 0, 0, 0, 0, t58 ^ 2 * t71, t61 * t58 * t82, t58 * t74, t61 ^ 2 * t71, t61 * t74, t63 ^ 2 * t87 (t48 * t63 - t55 * t79) * qJD(1) (-t49 * t63 + t55 * t81) * qJD(1) (-t48 * t58 + t49 * t61) * t77, t49 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1 + t55 ^ 2 / 0.2e1, t41 ^ 2 / 0.2e1, -t41 * t39, -t41 * t46, t39 ^ 2 / 0.2e1, t39 * t46, t46 ^ 2 / 0.2e1, -t20 * t46 + t26 * t39, t21 * t46 + t26 * t41, -t20 * t41 - t21 * t39, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t37, t29 ^ 2 / 0.2e1, -t29 * t37, t37 ^ 2 / 0.2e1, t18 * t29 + t9 * t37, -t10 * t37 + t18 * t31, -t10 * t29 - t9 * t31, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t22, -t84, t13, t75, -t83, t27, t7 * t23 + t3 * t28, t7 * t25 - t4 * t28, -t4 * t23 - t3 * t25, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t22, t13, t84, t27, t83, t75, -t1 * t28 + t5 * t23, t1 * t25 - t2 * t23, t2 * t28 - t5 * t25, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
