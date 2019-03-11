% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR9
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
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:55:42
% EndTime: 2019-03-09 22:55:42
% DurationCPUTime: 0.26s
% Computational Cost: add. (1716->72), mult. (3927->167), div. (0->0), fcn. (3136->12), ass. (0->56)
t64 = sin(qJ(2));
t65 = cos(qJ(2));
t59 = sin(pkin(6));
t74 = qJD(1) * t59;
t68 = t65 * t74;
t73 = cos(pkin(6)) * qJD(1);
t70 = pkin(1) * t73;
t46 = pkin(8) * t68 + t64 * t70;
t56 = qJD(2) + t73;
t39 = t56 * pkin(9) + t46;
t41 = (-pkin(2) * t65 - pkin(9) * t64 - pkin(1)) * t74;
t63 = sin(qJ(3));
t79 = cos(qJ(3));
t28 = -t63 * t39 + t79 * t41;
t69 = t64 * t74;
t44 = t63 * t56 + t79 * t69;
t51 = -qJD(3) + t68;
t20 = -t51 * pkin(3) - t44 * pkin(10) + t28;
t29 = t79 * t39 + t63 * t41;
t42 = -t79 * t56 + t63 * t69;
t23 = -t42 * pkin(10) + t29;
t62 = sin(qJ(4));
t78 = cos(qJ(4));
t12 = t62 * t20 + t78 * t23;
t48 = -qJD(4) + t51;
t10 = -t48 * qJ(5) + t12;
t31 = t78 * t42 + t62 * t44;
t33 = -t62 * t42 + t78 * t44;
t45 = -pkin(8) * t69 + t65 * t70;
t38 = -t56 * pkin(2) - t45;
t34 = t42 * pkin(3) + t38;
t15 = t31 * pkin(4) - t33 * qJ(5) + t34;
t58 = sin(pkin(12));
t75 = cos(pkin(12));
t6 = t75 * t10 + t58 * t15;
t77 = cos(qJ(6));
t66 = qJD(1) ^ 2;
t76 = t59 ^ 2 * t66;
t72 = t31 ^ 2 / 0.2e1;
t71 = t65 * t76;
t67 = t76 / 0.2e1;
t5 = -t58 * t10 + t75 * t15;
t11 = t78 * t20 - t62 * t23;
t9 = t48 * pkin(4) + qJD(5) - t11;
t61 = sin(qJ(6));
t30 = qJD(6) + t31;
t27 = t75 * t33 - t58 * t48;
t25 = t58 * t33 + t75 * t48;
t18 = -t61 * t25 + t77 * t27;
t16 = t77 * t25 + t61 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(11) + t6;
t3 = t31 * pkin(5) - t27 * pkin(11) + t5;
t2 = t61 * t3 + t77 * t4;
t1 = t77 * t3 - t61 * t4;
t8 = [0, 0, 0, 0, 0, t66 / 0.2e1, 0, 0, 0, 0, t64 ^ 2 * t67, t64 * t71, t56 * t69, t65 ^ 2 * t67, t56 * t68, t56 ^ 2 / 0.2e1, pkin(1) * t71 + t45 * t56, -pkin(1) * t64 * t76 - t46 * t56 (-t45 * t64 + t46 * t65) * t74, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t67, t44 ^ 2 / 0.2e1, -t44 * t42, -t44 * t51, t42 ^ 2 / 0.2e1, t42 * t51, t51 ^ 2 / 0.2e1, -t28 * t51 + t38 * t42, t29 * t51 + t38 * t44, -t28 * t44 - t29 * t42, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, -t33 * t48, t72, t31 * t48, t48 ^ 2 / 0.2e1, -t11 * t48 + t34 * t31, t12 * t48 + t34 * t33, -t11 * t33 - t12 * t31, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t31, t25 ^ 2 / 0.2e1, -t25 * t31, t72, t9 * t25 + t5 * t31, t9 * t27 - t6 * t31, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t30, t16 ^ 2 / 0.2e1, -t16 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t16, t7 * t18 - t2 * t30, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
