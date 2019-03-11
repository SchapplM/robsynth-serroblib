% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:43:08
% EndTime: 2019-03-10 01:43:08
% DurationCPUTime: 0.25s
% Computational Cost: add. (1295->68), mult. (3006->150), div. (0->0), fcn. (2356->10), ass. (0->55)
t71 = cos(pkin(6)) * qJD(1);
t56 = qJD(2) + t71;
t62 = sin(qJ(3));
t63 = sin(qJ(2));
t58 = sin(pkin(6));
t72 = qJD(1) * t58;
t68 = t63 * t72;
t76 = cos(qJ(3));
t42 = -t76 * t56 + t62 * t68;
t44 = t62 * t56 + t76 * t68;
t61 = sin(qJ(4));
t75 = cos(qJ(4));
t31 = t75 * t42 + t61 * t44;
t33 = -t61 * t42 + t75 * t44;
t64 = cos(qJ(2));
t69 = pkin(1) * t71;
t45 = -pkin(8) * t68 + t64 * t69;
t38 = -t56 * pkin(2) - t45;
t34 = t42 * pkin(3) + t38;
t13 = t31 * pkin(4) - t33 * pkin(11) + t34;
t60 = sin(qJ(5));
t74 = cos(qJ(5));
t67 = t64 * t72;
t46 = pkin(8) * t67 + t63 * t69;
t39 = t56 * pkin(9) + t46;
t41 = (-pkin(2) * t64 - pkin(9) * t63 - pkin(1)) * t72;
t27 = -t62 * t39 + t76 * t41;
t51 = -qJD(3) + t67;
t18 = -t51 * pkin(3) - t44 * pkin(10) + t27;
t28 = t76 * t39 + t62 * t41;
t21 = -t42 * pkin(10) + t28;
t10 = t61 * t18 + t75 * t21;
t48 = -qJD(4) + t51;
t8 = -t48 * pkin(11) + t10;
t4 = t60 * t13 + t74 * t8;
t65 = qJD(1) ^ 2;
t73 = t58 ^ 2 * t65;
t70 = t64 * t73;
t66 = t73 / 0.2e1;
t3 = t74 * t13 - t60 * t8;
t9 = t75 * t18 - t61 * t21;
t7 = t48 * pkin(4) - t9;
t30 = qJD(5) + t31;
t29 = t30 ^ 2 / 0.2e1;
t26 = t74 * t33 - t60 * t48;
t24 = t60 * t33 + t74 * t48;
t23 = t26 ^ 2 / 0.2e1;
t22 = t24 ^ 2 / 0.2e1;
t16 = t26 * t30;
t15 = t24 * t30;
t14 = t26 * t24;
t5 = t24 * pkin(5) + qJD(6) + t7;
t2 = -t24 * qJ(6) + t4;
t1 = t30 * pkin(5) - t26 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t65 / 0.2e1, 0, 0, 0, 0, t63 ^ 2 * t66, t63 * t70, t56 * t68, t64 ^ 2 * t66, t56 * t67, t56 ^ 2 / 0.2e1, pkin(1) * t70 + t45 * t56, -pkin(1) * t63 * t73 - t46 * t56 (-t45 * t63 + t46 * t64) * t72, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t66, t44 ^ 2 / 0.2e1, -t44 * t42, -t44 * t51, t42 ^ 2 / 0.2e1, t42 * t51, t51 ^ 2 / 0.2e1, -t27 * t51 + t38 * t42, t28 * t51 + t38 * t44, -t27 * t44 - t28 * t42, t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, -t33 * t48, t31 ^ 2 / 0.2e1, t31 * t48, t48 ^ 2 / 0.2e1, t34 * t31 - t9 * t48, t10 * t48 + t34 * t33, -t10 * t31 - t9 * t33, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t23, -t14, t16, t22, -t15, t29, t7 * t24 + t3 * t30, t7 * t26 - t4 * t30, -t4 * t24 - t3 * t26, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t23, -t14, t16, t22, -t15, t29, t1 * t30 + t5 * t24, -t2 * t30 + t5 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
