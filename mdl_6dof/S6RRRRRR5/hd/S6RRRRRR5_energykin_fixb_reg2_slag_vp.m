% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:04:46
% EndTime: 2019-03-10 04:04:46
% DurationCPUTime: 0.26s
% Computational Cost: add. (1809->72), mult. (4223->170), div. (0->0), fcn. (3384->12), ass. (0->56)
t65 = sin(qJ(2));
t66 = cos(qJ(2));
t59 = sin(pkin(6));
t74 = qJD(1) * t59;
t69 = t66 * t74;
t73 = cos(pkin(6)) * qJD(1);
t71 = pkin(1) * t73;
t45 = pkin(8) * t69 + t65 * t71;
t57 = qJD(2) + t73;
t38 = t57 * pkin(9) + t45;
t40 = (-pkin(2) * t66 - pkin(9) * t65 - pkin(1)) * t74;
t64 = sin(qJ(3));
t79 = cos(qJ(3));
t27 = -t64 * t38 + t79 * t40;
t70 = t65 * t74;
t43 = t64 * t57 + t79 * t70;
t52 = -qJD(3) + t69;
t24 = -t52 * pkin(3) - t43 * pkin(10) + t27;
t28 = t79 * t38 + t64 * t40;
t41 = -t79 * t57 + t64 * t70;
t26 = -t41 * pkin(10) + t28;
t63 = sin(qJ(4));
t78 = cos(qJ(4));
t13 = t63 * t24 + t78 * t26;
t30 = t78 * t41 + t63 * t43;
t11 = -t30 * pkin(11) + t13;
t62 = sin(qJ(5));
t77 = cos(qJ(5));
t12 = t78 * t24 - t63 * t26;
t32 = -t63 * t41 + t78 * t43;
t49 = -qJD(4) + t52;
t9 = -t49 * pkin(4) - t32 * pkin(11) + t12;
t6 = t77 * t11 + t62 * t9;
t76 = cos(qJ(6));
t67 = qJD(1) ^ 2;
t75 = t59 ^ 2 * t67;
t72 = t66 * t75;
t68 = t75 / 0.2e1;
t18 = t77 * t30 + t62 * t32;
t44 = -pkin(8) * t70 + t66 * t71;
t5 = -t62 * t11 + t77 * t9;
t37 = -t57 * pkin(2) - t44;
t33 = t41 * pkin(3) + t37;
t21 = t30 * pkin(4) + t33;
t61 = sin(qJ(6));
t47 = -qJD(5) + t49;
t20 = -t62 * t30 + t77 * t32;
t17 = qJD(6) + t18;
t16 = t76 * t20 - t61 * t47;
t14 = t61 * t20 + t76 * t47;
t7 = t18 * pkin(5) - t20 * pkin(12) + t21;
t4 = -t47 * pkin(12) + t6;
t3 = t47 * pkin(5) - t5;
t2 = t76 * t4 + t61 * t7;
t1 = -t61 * t4 + t76 * t7;
t8 = [0, 0, 0, 0, 0, t67 / 0.2e1, 0, 0, 0, 0, t65 ^ 2 * t68, t65 * t72, t57 * t70, t66 ^ 2 * t68, t57 * t69, t57 ^ 2 / 0.2e1, pkin(1) * t72 + t44 * t57, -pkin(1) * t65 * t75 - t45 * t57 (-t44 * t65 + t45 * t66) * t74, t45 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t68, t43 ^ 2 / 0.2e1, -t43 * t41, -t43 * t52, t41 ^ 2 / 0.2e1, t41 * t52, t52 ^ 2 / 0.2e1, -t27 * t52 + t37 * t41, t28 * t52 + t37 * t43, -t27 * t43 - t28 * t41, t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t30, -t32 * t49, t30 ^ 2 / 0.2e1, t30 * t49, t49 ^ 2 / 0.2e1, -t12 * t49 + t33 * t30, t13 * t49 + t33 * t32, -t12 * t32 - t13 * t30, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, -t20 * t47, t18 ^ 2 / 0.2e1, t18 * t47, t47 ^ 2 / 0.2e1, t21 * t18 - t5 * t47, t21 * t20 + t6 * t47, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t17, t14 ^ 2 / 0.2e1, -t14 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t14, t3 * t16 - t2 * t17, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
