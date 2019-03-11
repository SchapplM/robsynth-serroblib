% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:52:55
% EndTime: 2019-03-09 05:52:56
% DurationCPUTime: 0.31s
% Computational Cost: add. (1454->79), mult. (4731->177), div. (0->0), fcn. (3877->12), ass. (0->68)
t55 = sin(pkin(12));
t58 = cos(pkin(12));
t57 = sin(pkin(6));
t78 = qJD(1) * t57;
t70 = qJ(2) * t78;
t60 = cos(pkin(6));
t77 = qJD(1) * t60;
t73 = pkin(1) * t77;
t46 = t55 * t73 + t58 * t70;
t56 = sin(pkin(7));
t59 = cos(pkin(7));
t80 = t57 * t58;
t32 = (t56 * t60 + t59 * t80) * qJD(1) * pkin(9) + t46;
t51 = t58 * t73;
t82 = t55 * t57;
t35 = t51 + (pkin(2) * t60 + (-pkin(9) * t59 - qJ(2)) * t82) * qJD(1);
t41 = qJD(2) + (-pkin(9) * t55 * t56 - pkin(2) * t58 - pkin(1)) * t78;
t63 = sin(qJ(3));
t65 = cos(qJ(3));
t17 = -t63 * t32 + (t35 * t59 + t41 * t56) * t65;
t66 = qJD(1) ^ 2;
t89 = t66 / 0.2e1;
t88 = pkin(4) + pkin(11);
t87 = cos(qJ(6));
t79 = t59 * t63;
t81 = t56 * t63;
t38 = (t60 * t81 + (t55 * t65 + t58 * t79) * t57) * qJD(1);
t72 = t58 * t78;
t43 = t56 * t72 - t59 * t77 - qJD(3);
t62 = sin(qJ(4));
t64 = cos(qJ(4));
t25 = t62 * t38 + t64 * t43;
t27 = t64 * t38 - t62 * t43;
t86 = t27 * t25;
t36 = t63 * t55 * t78 + (-t56 * t77 - t59 * t72) * t65;
t34 = qJD(4) + t36;
t85 = t27 * t34;
t84 = t34 * t25;
t83 = t57 ^ 2 * t66;
t22 = -t56 * t35 + t59 * t41;
t12 = t36 * pkin(3) - t38 * pkin(10) + t22;
t18 = t65 * t32 + t35 * t79 + t41 * t81;
t16 = -t43 * pkin(10) + t18;
t9 = t62 * t12 + t64 * t16;
t76 = t25 ^ 2 / 0.2e1;
t75 = t27 ^ 2 / 0.2e1;
t74 = t57 * t60 * t66;
t71 = t83 / 0.2e1;
t8 = t64 * t12 - t62 * t16;
t7 = -t34 * qJ(5) - t9;
t69 = qJD(5) - t8;
t15 = t43 * pkin(3) - t17;
t67 = -t27 * qJ(5) + t15;
t61 = sin(qJ(6));
t52 = -pkin(1) * t78 + qJD(2);
t45 = -t55 * t70 + t51;
t33 = t34 ^ 2 / 0.2e1;
t24 = qJD(6) + t27;
t21 = t61 * t25 + t87 * t34;
t19 = -t87 * t25 + t61 * t34;
t10 = t25 * pkin(4) + t67;
t6 = -t34 * pkin(4) + t69;
t5 = t88 * t25 + t67;
t4 = -t25 * pkin(5) - t7;
t3 = t27 * pkin(5) - t88 * t34 + t69;
t2 = t61 * t3 + t87 * t5;
t1 = t87 * t3 - t61 * t5;
t11 = [0, 0, 0, 0, 0, t89, 0, 0, 0, 0, t55 ^ 2 * t71, t55 * t58 * t83, t55 * t74, t58 ^ 2 * t71, t58 * t74, t60 ^ 2 * t89 (t45 * t60 - t52 * t80) * qJD(1) (-t46 * t60 + t52 * t82) * qJD(1) (-t45 * t55 + t46 * t58) * t78, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1, t38 ^ 2 / 0.2e1, -t38 * t36, -t38 * t43, t36 ^ 2 / 0.2e1, t36 * t43, t43 ^ 2 / 0.2e1, -t17 * t43 + t22 * t36, t18 * t43 + t22 * t38, -t17 * t38 - t18 * t36, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t75, -t86, t85, t76, -t84, t33, t15 * t25 + t8 * t34, t15 * t27 - t9 * t34, -t9 * t25 - t8 * t27, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t33, -t85, t84, t75, -t86, t76, t7 * t25 + t6 * t27, -t10 * t25 + t6 * t34, -t10 * t27 - t7 * t34, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * t24, t19 ^ 2 / 0.2e1, -t19 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t4 * t19, -t2 * t24 + t4 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
