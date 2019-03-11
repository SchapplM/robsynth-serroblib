% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:23:32
% EndTime: 2019-03-09 23:23:33
% DurationCPUTime: 0.26s
% Computational Cost: add. (1790->72), mult. (4079->167), div. (0->0), fcn. (3255->12), ass. (0->56)
t66 = cos(qJ(2));
t64 = sin(qJ(2));
t59 = sin(pkin(6));
t74 = qJD(1) * t59;
t70 = t64 * t74;
t73 = cos(pkin(6)) * qJD(1);
t71 = pkin(1) * t73;
t47 = -pkin(8) * t70 + t66 * t71;
t56 = qJD(2) + t73;
t38 = -pkin(2) * t56 - t47;
t63 = sin(qJ(3));
t65 = cos(qJ(3));
t44 = -t56 * t65 + t63 * t70;
t46 = t56 * t63 + t65 * t70;
t25 = pkin(3) * t44 - pkin(10) * t46 + t38;
t69 = t66 * t74;
t48 = pkin(8) * t69 + t64 * t71;
t39 = pkin(9) * t56 + t48;
t42 = (-pkin(2) * t66 - pkin(9) * t64 - pkin(1)) * t74;
t30 = t39 * t65 + t42 * t63;
t51 = -qJD(3) + t69;
t28 = -pkin(10) * t51 + t30;
t62 = sin(qJ(4));
t78 = cos(qJ(4));
t17 = t25 * t62 + t28 * t78;
t32 = t46 * t62 + t51 * t78;
t12 = -qJ(5) * t32 + t17;
t58 = sin(pkin(12));
t75 = cos(pkin(12));
t16 = t25 * t78 - t28 * t62;
t34 = t46 * t78 - t51 * t62;
t43 = qJD(4) + t44;
t9 = pkin(4) * t43 - qJ(5) * t34 + t16;
t6 = t12 * t75 + t58 * t9;
t77 = cos(qJ(6));
t67 = qJD(1) ^ 2;
t76 = t59 ^ 2 * t67;
t72 = t66 * t76;
t68 = t76 / 0.2e1;
t5 = -t12 * t58 + t75 * t9;
t29 = -t39 * t63 + t42 * t65;
t27 = pkin(3) * t51 - t29;
t18 = pkin(4) * t32 + qJD(5) + t27;
t61 = sin(qJ(6));
t41 = qJD(6) + t43;
t40 = t43 ^ 2 / 0.2e1;
t22 = -t32 * t58 + t34 * t75;
t20 = t32 * t75 + t34 * t58;
t15 = -t20 * t61 + t22 * t77;
t13 = t20 * t77 + t22 * t61;
t10 = pkin(5) * t20 + t18;
t4 = -pkin(11) * t20 + t6;
t3 = pkin(5) * t43 - pkin(11) * t22 + t5;
t2 = t3 * t61 + t4 * t77;
t1 = t3 * t77 - t4 * t61;
t7 = [0, 0, 0, 0, 0, t67 / 0.2e1, 0, 0, 0, 0, t64 ^ 2 * t68, t64 * t72, t56 * t70, t66 ^ 2 * t68, t56 * t69, t56 ^ 2 / 0.2e1, pkin(1) * t72 + t47 * t56, -pkin(1) * t64 * t76 - t48 * t56 (-t47 * t64 + t48 * t66) * t74, t48 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t68, t46 ^ 2 / 0.2e1, -t46 * t44, -t46 * t51, t44 ^ 2 / 0.2e1, t44 * t51, t51 ^ 2 / 0.2e1, -t29 * t51 + t38 * t44, t30 * t51 + t38 * t46, -t29 * t46 - t30 * t44, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t43, t32 ^ 2 / 0.2e1, -t32 * t43, t40, t16 * t43 + t27 * t32, -t17 * t43 + t27 * t34, -t16 * t34 - t17 * t32, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t43, t20 ^ 2 / 0.2e1, -t20 * t43, t40, t18 * t20 + t43 * t5, t18 * t22 - t43 * t6, -t20 * t6 - t22 * t5, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t41, t13 ^ 2 / 0.2e1, -t13 * t41, t41 ^ 2 / 0.2e1, t1 * t41 + t10 * t13, t10 * t15 - t2 * t41, -t1 * t15 - t13 * t2, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t7;
