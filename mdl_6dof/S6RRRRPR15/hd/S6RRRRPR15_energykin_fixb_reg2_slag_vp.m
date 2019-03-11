% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR15_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:50:11
% EndTime: 2019-03-10 00:50:11
% DurationCPUTime: 0.33s
% Computational Cost: add. (1762->77), mult. (4732->166), div. (0->0), fcn. (3875->12), ass. (0->65)
t78 = cos(pkin(6)) * qJD(1);
t53 = qJD(2) + t78;
t55 = sin(pkin(7));
t57 = cos(pkin(7));
t65 = cos(qJ(2));
t56 = sin(pkin(6));
t79 = qJD(1) * t56;
t72 = t65 * t79;
t89 = t53 * t55 + t57 * t72;
t62 = sin(qJ(2));
t74 = pkin(1) * t78;
t46 = pkin(9) * t72 + t62 * t74;
t33 = t89 * pkin(10) + t46;
t52 = t65 * t74;
t73 = t62 * t79;
t35 = t53 * pkin(2) + t52 + (-pkin(10) * t57 - pkin(9)) * t73;
t42 = (-pkin(10) * t55 * t62 - pkin(2) * t65 - pkin(1)) * t79;
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t17 = -t61 * t33 + (t35 * t57 + t42 * t55) * t64;
t88 = pkin(4) + pkin(12);
t87 = cos(qJ(6));
t80 = t57 * t61;
t81 = t55 * t61;
t38 = t53 * t81 + (t62 * t64 + t65 * t80) * t79;
t43 = -t57 * t53 + t55 * t72 - qJD(3);
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t25 = t60 * t38 + t63 * t43;
t27 = t63 * t38 - t60 * t43;
t86 = t27 * t25;
t36 = t61 * t73 - t89 * t64;
t34 = qJD(4) + t36;
t85 = t27 * t34;
t84 = t34 * t25;
t66 = qJD(1) ^ 2;
t82 = t56 ^ 2 * t66;
t22 = -t55 * t35 + t57 * t42;
t12 = t36 * pkin(3) - t38 * pkin(11) + t22;
t18 = t64 * t33 + t35 * t80 + t42 * t81;
t16 = -t43 * pkin(11) + t18;
t9 = t60 * t12 + t63 * t16;
t77 = t25 ^ 2 / 0.2e1;
t76 = t27 ^ 2 / 0.2e1;
t75 = t65 * t82;
t71 = t82 / 0.2e1;
t8 = t63 * t12 - t60 * t16;
t7 = -t34 * qJ(5) - t9;
t69 = qJD(5) - t8;
t15 = t43 * pkin(3) - t17;
t67 = -t27 * qJ(5) + t15;
t59 = sin(qJ(6));
t45 = -pkin(9) * t73 + t52;
t32 = t34 ^ 2 / 0.2e1;
t24 = qJD(6) + t27;
t21 = t59 * t25 + t87 * t34;
t19 = -t87 * t25 + t59 * t34;
t10 = t25 * pkin(4) + t67;
t6 = -t34 * pkin(4) + t69;
t5 = t88 * t25 + t67;
t4 = -t25 * pkin(5) - t7;
t3 = t27 * pkin(5) - t88 * t34 + t69;
t2 = t59 * t3 + t87 * t5;
t1 = t87 * t3 - t59 * t5;
t11 = [0, 0, 0, 0, 0, t66 / 0.2e1, 0, 0, 0, 0, t62 ^ 2 * t71, t62 * t75, t53 * t73, t65 ^ 2 * t71, t53 * t72, t53 ^ 2 / 0.2e1, pkin(1) * t75 + t45 * t53, -pkin(1) * t62 * t82 - t46 * t53 (-t45 * t62 + t46 * t65) * t79, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t71, t38 ^ 2 / 0.2e1, -t38 * t36, -t38 * t43, t36 ^ 2 / 0.2e1, t36 * t43, t43 ^ 2 / 0.2e1, -t17 * t43 + t22 * t36, t18 * t43 + t22 * t38, -t17 * t38 - t18 * t36, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t76, -t86, t85, t77, -t84, t32, t15 * t25 + t8 * t34, t15 * t27 - t9 * t34, -t9 * t25 - t8 * t27, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t32, -t85, t84, t76, -t86, t77, t7 * t25 + t6 * t27, -t10 * t25 + t6 * t34, -t10 * t27 - t7 * t34, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * t24, t19 ^ 2 / 0.2e1, -t19 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t4 * t19, -t2 * t24 + t4 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
