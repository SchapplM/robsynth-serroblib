% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_energykin_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:55:52
% EndTime: 2019-03-08 20:55:52
% DurationCPUTime: 0.32s
% Computational Cost: add. (1348->72), mult. (4109->185), div. (0->0), fcn. (3456->16), ass. (0->66)
t68 = cos(qJ(2));
t58 = sin(pkin(6));
t81 = qJD(1) * t58;
t47 = qJD(2) * pkin(2) + t68 * t81;
t57 = sin(pkin(7));
t61 = cos(pkin(7));
t62 = cos(pkin(6));
t80 = qJD(1) * t62;
t92 = t47 * t61 + t57 * t80;
t66 = sin(qJ(2));
t79 = qJD(2) * t57;
t45 = qJ(3) * t79 + t66 * t81;
t55 = sin(pkin(14));
t59 = cos(pkin(14));
t31 = t59 * t45 + t92 * t55;
t56 = sin(pkin(8));
t60 = cos(pkin(8));
t83 = t57 * t60;
t22 = (t56 * t61 + t59 * t83) * qJD(2) * pkin(10) + t31;
t30 = -t55 * t45 + t92 * t59;
t90 = pkin(10) * t55;
t23 = (pkin(3) * t61 - t83 * t90) * qJD(2) + t30;
t77 = t61 * t80 + qJD(3);
t32 = (-t47 + (-pkin(3) * t59 - t56 * t90) * qJD(2)) * t57 + t77;
t65 = sin(qJ(4));
t67 = cos(qJ(4));
t12 = -t65 * t22 + (t23 * t60 + t32 * t56) * t67;
t69 = qJD(2) ^ 2;
t91 = t69 / 0.2e1;
t82 = t60 * t65;
t84 = t56 * t65;
t13 = t67 * t22 + t23 * t82 + t32 * t84;
t74 = t59 * t79;
t78 = qJD(2) * t61;
t40 = t56 * t74 - t60 * t78 - qJD(4);
t10 = -t40 * pkin(11) + t13;
t15 = -t56 * t23 + t60 * t32;
t35 = t65 * t55 * t79 + (-t56 * t78 - t60 * t74) * t67;
t37 = (t61 * t84 + (t55 * t67 + t59 * t82) * t57) * qJD(2);
t14 = t35 * pkin(4) - t37 * pkin(11) + t15;
t64 = sin(qJ(5));
t89 = cos(qJ(5));
t6 = t89 * t10 + t64 * t14;
t88 = cos(qJ(6));
t38 = -t57 * t47 + t77;
t87 = t38 * t57;
t85 = t57 ^ 2 * t69;
t76 = t57 * t61 * t69;
t73 = t85 / 0.2e1;
t72 = qJD(2) * t81;
t25 = t64 * t37 + t89 * t40;
t5 = -t64 * t10 + t89 * t14;
t9 = t40 * pkin(4) - t12;
t70 = qJD(1) ^ 2;
t63 = sin(qJ(6));
t34 = qJD(5) + t35;
t27 = t89 * t37 - t64 * t40;
t24 = qJD(6) + t25;
t18 = t88 * t27 + t63 * t34;
t16 = t63 * t27 - t88 * t34;
t7 = t25 * pkin(5) - t27 * pkin(12) + t9;
t4 = t34 * pkin(12) + t6;
t3 = -t34 * pkin(5) - t5;
t2 = t88 * t4 + t63 * t7;
t1 = -t63 * t4 + t88 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t70 / 0.2e1, 0, 0, 0, 0, 0, t91, t68 * t72, -t66 * t72, 0 (t62 ^ 2 / 0.2e1 + (t66 ^ 2 / 0.2e1 + t68 ^ 2 / 0.2e1) * t58 ^ 2) * t70, t55 ^ 2 * t73, t55 * t59 * t85, t55 * t76, t59 ^ 2 * t73, t59 * t76, t61 ^ 2 * t91 (t30 * t61 - t59 * t87) * qJD(2) (-t31 * t61 + t55 * t87) * qJD(2) (-t30 * t55 + t31 * t59) * t79, t31 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t37 ^ 2 / 0.2e1, -t35 * t37, -t40 * t37, t35 ^ 2 / 0.2e1, t40 * t35, t40 ^ 2 / 0.2e1, -t12 * t40 + t15 * t35, t13 * t40 + t15 * t37, -t12 * t37 - t13 * t35, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t34, t25 ^ 2 / 0.2e1, -t25 * t34, t34 ^ 2 / 0.2e1, t9 * t25 + t5 * t34, t9 * t27 - t6 * t34, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t24, t16 ^ 2 / 0.2e1, -t16 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t3 * t16, t3 * t18 - t2 * t24, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
