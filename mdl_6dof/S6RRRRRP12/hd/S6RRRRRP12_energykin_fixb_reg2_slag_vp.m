% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:20:36
% EndTime: 2019-03-10 03:20:36
% DurationCPUTime: 0.34s
% Computational Cost: add. (2021->74), mult. (5403->166), div. (0->0), fcn. (4447->12), ass. (0->62)
t77 = cos(pkin(6)) * qJD(1);
t56 = qJD(2) + t77;
t58 = sin(pkin(7));
t60 = cos(pkin(7));
t67 = cos(qJ(2));
t59 = sin(pkin(6));
t78 = qJD(1) * t59;
t72 = t67 * t78;
t87 = t56 * t58 + t60 * t72;
t65 = sin(qJ(2));
t74 = pkin(1) * t77;
t49 = pkin(9) * t72 + t65 * t74;
t36 = t87 * pkin(10) + t49;
t55 = t67 * t74;
t73 = t65 * t78;
t38 = t56 * pkin(2) + t55 + (-pkin(10) * t60 - pkin(9)) * t73;
t45 = (-pkin(10) * t58 * t65 - pkin(2) * t67 - pkin(1)) * t78;
t64 = sin(qJ(3));
t66 = cos(qJ(3));
t20 = -t64 * t36 + (t38 * t60 + t45 * t58) * t66;
t46 = -t60 * t56 + t58 * t72 - qJD(3);
t18 = t46 * pkin(3) - t20;
t79 = t60 * t64;
t80 = t58 * t64;
t41 = t56 * t80 + (t65 * t66 + t67 * t79) * t78;
t63 = sin(qJ(4));
t86 = cos(qJ(4));
t29 = t63 * t41 + t86 * t46;
t31 = t41 * t86 - t63 * t46;
t12 = t29 * pkin(4) - t31 * pkin(12) + t18;
t62 = sin(qJ(5));
t26 = -t58 * t38 + t60 * t45;
t39 = t64 * t73 - t87 * t66;
t15 = t39 * pkin(3) - t41 * pkin(11) + t26;
t21 = t66 * t36 + t38 * t79 + t45 * t80;
t19 = -t46 * pkin(11) + t21;
t10 = t63 * t15 + t86 * t19;
t37 = qJD(4) + t39;
t8 = t37 * pkin(12) + t10;
t85 = cos(qJ(5));
t4 = t62 * t12 + t85 * t8;
t23 = t62 * t31 - t37 * t85;
t25 = t31 * t85 + t62 * t37;
t84 = t25 * t23;
t28 = qJD(5) + t29;
t83 = t28 * t23;
t68 = qJD(1) ^ 2;
t81 = t59 ^ 2 * t68;
t76 = t23 ^ 2 / 0.2e1;
t75 = t67 * t81;
t71 = t81 / 0.2e1;
t9 = t15 * t86 - t63 * t19;
t3 = t12 * t85 - t62 * t8;
t7 = -t37 * pkin(4) - t9;
t48 = -pkin(9) * t73 + t55;
t27 = t28 ^ 2 / 0.2e1;
t22 = t25 ^ 2 / 0.2e1;
t13 = t25 * t28;
t5 = t23 * pkin(5) - t25 * qJ(6) + t7;
t2 = t28 * qJ(6) + t4;
t1 = -t28 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t68 / 0.2e1, 0, 0, 0, 0, t65 ^ 2 * t71, t65 * t75, t56 * t73, t67 ^ 2 * t71, t56 * t72, t56 ^ 2 / 0.2e1, pkin(1) * t75 + t48 * t56, -pkin(1) * t65 * t81 - t49 * t56 (-t48 * t65 + t49 * t67) * t78, t49 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t71, t41 ^ 2 / 0.2e1, -t41 * t39, -t41 * t46, t39 ^ 2 / 0.2e1, t39 * t46, t46 ^ 2 / 0.2e1, -t20 * t46 + t26 * t39, t21 * t46 + t26 * t41, -t20 * t41 - t21 * t39, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t37, t29 ^ 2 / 0.2e1, -t29 * t37, t37 ^ 2 / 0.2e1, t18 * t29 + t9 * t37, -t10 * t37 + t18 * t31, -t10 * t29 - t9 * t31, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t22, -t84, t13, t76, -t83, t27, t7 * t23 + t3 * t28, t7 * t25 - t4 * t28, -t4 * t23 - t3 * t25, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t22, t13, t84, t27, t83, t76, -t1 * t28 + t5 * t23, t1 * t25 - t2 * t23, t2 * t28 - t5 * t25, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
