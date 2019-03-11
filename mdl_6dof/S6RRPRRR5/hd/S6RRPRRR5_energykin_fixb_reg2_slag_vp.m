% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:46:46
% EndTime: 2019-03-09 13:46:46
% DurationCPUTime: 0.26s
% Computational Cost: add. (1418->73), mult. (3892->168), div. (0->0), fcn. (3104->12), ass. (0->57)
t66 = cos(qJ(2));
t73 = cos(pkin(6)) * qJD(1);
t71 = pkin(1) * t73;
t54 = t66 * t71;
t55 = qJD(2) + t73;
t64 = sin(qJ(2));
t58 = sin(pkin(6));
t74 = qJD(1) * t58;
t70 = t64 * t74;
t37 = t55 * pkin(2) + t54 + (-pkin(8) - qJ(3)) * t70;
t69 = t66 * t74;
t47 = pkin(8) * t69 + t64 * t71;
t40 = qJ(3) * t69 + t47;
t57 = sin(pkin(12));
t59 = cos(pkin(12));
t25 = t57 * t37 + t59 * t40;
t22 = t55 * pkin(9) + t25;
t43 = t57 * t70 - t59 * t69;
t45 = (t57 * t66 + t59 * t64) * t74;
t48 = qJD(3) + (-pkin(2) * t66 - pkin(1)) * t74;
t30 = t43 * pkin(3) - t45 * pkin(9) + t48;
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t15 = t65 * t22 + t63 * t30;
t42 = qJD(4) + t43;
t10 = t42 * pkin(10) + t15;
t24 = t59 * t37 - t57 * t40;
t21 = -t55 * pkin(3) - t24;
t34 = t63 * t45 - t65 * t55;
t36 = t65 * t45 + t63 * t55;
t13 = t34 * pkin(4) - t36 * pkin(10) + t21;
t62 = sin(qJ(5));
t77 = cos(qJ(5));
t6 = t77 * t10 + t62 * t13;
t76 = cos(qJ(6));
t67 = qJD(1) ^ 2;
t75 = t58 ^ 2 * t67;
t72 = t66 * t75;
t68 = t75 / 0.2e1;
t5 = -t62 * t10 + t77 * t13;
t14 = -t63 * t22 + t65 * t30;
t33 = qJD(5) + t34;
t9 = -t42 * pkin(4) - t14;
t61 = sin(qJ(6));
t51 = t55 ^ 2 / 0.2e1;
t46 = -pkin(8) * t70 + t54;
t31 = qJD(6) + t33;
t28 = t77 * t36 + t62 * t42;
t26 = t62 * t36 - t77 * t42;
t18 = -t61 * t26 + t76 * t28;
t16 = t76 * t26 + t61 * t28;
t7 = t26 * pkin(5) + t9;
t4 = -t26 * pkin(11) + t6;
t3 = t33 * pkin(5) - t28 * pkin(11) + t5;
t2 = t61 * t3 + t76 * t4;
t1 = t76 * t3 - t61 * t4;
t8 = [0, 0, 0, 0, 0, t67 / 0.2e1, 0, 0, 0, 0, t64 ^ 2 * t68, t64 * t72, t55 * t70, t66 ^ 2 * t68, t55 * t69, t51, pkin(1) * t72 + t46 * t55, -pkin(1) * t64 * t75 - t47 * t55 (-t46 * t64 + t47 * t66) * t74, t47 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t68, t45 ^ 2 / 0.2e1, -t45 * t43, t45 * t55, t43 ^ 2 / 0.2e1, -t43 * t55, t51, t24 * t55 + t48 * t43, -t25 * t55 + t48 * t45, -t24 * t45 - t25 * t43, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1, t36 ^ 2 / 0.2e1, -t36 * t34, t36 * t42, t34 ^ 2 / 0.2e1, -t34 * t42, t42 ^ 2 / 0.2e1, t14 * t42 + t21 * t34, -t15 * t42 + t21 * t36, -t14 * t36 - t15 * t34, t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * t33, t26 ^ 2 / 0.2e1, -t26 * t33, t33 ^ 2 / 0.2e1, t9 * t26 + t5 * t33, t9 * t28 - t6 * t33, -t6 * t26 - t5 * t28, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t31, t16 ^ 2 / 0.2e1, -t16 * t31, t31 ^ 2 / 0.2e1, t1 * t31 + t7 * t16, t7 * t18 - t2 * t31, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
