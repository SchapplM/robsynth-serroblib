% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRP6
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:01:07
% EndTime: 2019-03-09 17:01:07
% DurationCPUTime: 0.23s
% Computational Cost: add. (1258->68), mult. (3006->147), div. (0->0), fcn. (2356->10), ass. (0->55)
t63 = cos(qJ(2));
t62 = sin(qJ(2));
t58 = sin(pkin(6));
t71 = qJD(1) * t58;
t67 = t62 * t71;
t70 = cos(pkin(6)) * qJD(1);
t68 = pkin(1) * t70;
t45 = -pkin(8) * t67 + t63 * t68;
t55 = qJD(2) + t70;
t38 = -t55 * pkin(2) - t45;
t61 = sin(qJ(3));
t75 = cos(qJ(3));
t42 = -t75 * t55 + t61 * t67;
t31 = t42 * pkin(3) + qJD(4) + t38;
t44 = t61 * t55 + t75 * t67;
t57 = sin(pkin(11));
t72 = cos(pkin(11));
t32 = t72 * t42 + t57 * t44;
t34 = -t57 * t42 + t72 * t44;
t13 = t32 * pkin(4) - t34 * pkin(10) + t31;
t60 = sin(qJ(5));
t74 = cos(qJ(5));
t66 = t63 * t71;
t46 = pkin(8) * t66 + t62 * t68;
t39 = t55 * pkin(9) + t46;
t41 = (-pkin(2) * t63 - pkin(9) * t62 - pkin(1)) * t71;
t27 = -t61 * t39 + t75 * t41;
t50 = -qJD(3) + t66;
t18 = -t50 * pkin(3) - t44 * qJ(4) + t27;
t28 = t75 * t39 + t61 * t41;
t21 = -t42 * qJ(4) + t28;
t10 = t57 * t18 + t72 * t21;
t8 = -t50 * pkin(10) + t10;
t4 = t60 * t13 + t74 * t8;
t64 = qJD(1) ^ 2;
t73 = t58 ^ 2 * t64;
t69 = t63 * t73;
t65 = t73 / 0.2e1;
t3 = t74 * t13 - t60 * t8;
t9 = t72 * t18 - t57 * t21;
t7 = t50 * pkin(4) - t9;
t48 = t50 ^ 2 / 0.2e1;
t30 = qJD(5) + t32;
t29 = t30 ^ 2 / 0.2e1;
t26 = t74 * t34 - t60 * t50;
t24 = t60 * t34 + t74 * t50;
t23 = t26 ^ 2 / 0.2e1;
t22 = t24 ^ 2 / 0.2e1;
t16 = t26 * t30;
t15 = t24 * t30;
t14 = t26 * t24;
t5 = t24 * pkin(5) + qJD(6) + t7;
t2 = -t24 * qJ(6) + t4;
t1 = t30 * pkin(5) - t26 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t64 / 0.2e1, 0, 0, 0, 0, t62 ^ 2 * t65, t62 * t69, t55 * t67, t63 ^ 2 * t65, t55 * t66, t55 ^ 2 / 0.2e1, pkin(1) * t69 + t45 * t55, -pkin(1) * t62 * t73 - t46 * t55 (-t45 * t62 + t46 * t63) * t71, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t65, t44 ^ 2 / 0.2e1, -t44 * t42, -t44 * t50, t42 ^ 2 / 0.2e1, t42 * t50, t48, -t27 * t50 + t38 * t42, t28 * t50 + t38 * t44, -t27 * t44 - t28 * t42, t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, -t34 * t50, t32 ^ 2 / 0.2e1, t32 * t50, t48, t31 * t32 - t9 * t50, t10 * t50 + t31 * t34, -t10 * t32 - t9 * t34, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t23, -t14, t16, t22, -t15, t29, t7 * t24 + t3 * t30, t7 * t26 - t4 * t30, -t4 * t24 - t3 * t26, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t23, -t14, t16, t22, -t15, t29, t1 * t30 + t5 * t24, -t2 * t30 + t5 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
