% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:39:11
% EndTime: 2019-03-10 05:39:12
% DurationCPUTime: 0.41s
% Computational Cost: add. (2696->80), mult. (7111->186), div. (0->0), fcn. (5930->14), ass. (0->63)
t83 = cos(pkin(6)) * qJD(1);
t61 = qJD(2) + t83;
t63 = sin(pkin(7));
t65 = cos(pkin(7));
t74 = cos(qJ(2));
t64 = sin(pkin(6));
t84 = qJD(1) * t64;
t79 = t74 * t84;
t91 = t61 * t63 + t65 * t79;
t71 = sin(qJ(2));
t81 = pkin(1) * t83;
t54 = pkin(9) * t79 + t71 * t81;
t41 = t91 * pkin(10) + t54;
t60 = t74 * t81;
t80 = t71 * t84;
t43 = t61 * pkin(2) + t60 + (-pkin(10) * t65 - pkin(9)) * t80;
t50 = (-pkin(10) * t63 * t71 - pkin(2) * t74 - pkin(1)) * t84;
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t25 = -t70 * t41 + (t43 * t65 + t50 * t63) * t73;
t31 = -t63 * t43 + t65 * t50;
t44 = t70 * t80 - t91 * t73;
t85 = t65 * t70;
t86 = t63 * t70;
t46 = t61 * t86 + (t71 * t73 + t74 * t85) * t84;
t20 = t44 * pkin(3) - t46 * pkin(11) + t31;
t26 = t73 * t41 + t43 * t85 + t50 * t86;
t51 = -t65 * t61 + t63 * t79 - qJD(3);
t24 = -t51 * pkin(11) + t26;
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t12 = t69 * t20 + t72 * t24;
t42 = qJD(4) + t44;
t10 = t42 * pkin(12) + t12;
t23 = t51 * pkin(3) - t25;
t34 = t69 * t46 + t72 * t51;
t36 = t72 * t46 - t69 * t51;
t15 = t34 * pkin(4) - t36 * pkin(12) + t23;
t68 = sin(qJ(5));
t90 = cos(qJ(5));
t6 = t90 * t10 + t68 * t15;
t89 = cos(qJ(6));
t75 = qJD(1) ^ 2;
t87 = t64 ^ 2 * t75;
t82 = t74 * t87;
t78 = t87 / 0.2e1;
t5 = -t68 * t10 + t90 * t15;
t11 = t72 * t20 - t69 * t24;
t33 = qJD(5) + t34;
t9 = -t42 * pkin(4) - t11;
t67 = sin(qJ(6));
t53 = -pkin(9) * t80 + t60;
t32 = qJD(6) + t33;
t30 = t90 * t36 + t68 * t42;
t28 = t68 * t36 - t90 * t42;
t18 = -t67 * t28 + t89 * t30;
t16 = t89 * t28 + t67 * t30;
t7 = t28 * pkin(5) + t9;
t4 = -t28 * pkin(13) + t6;
t3 = t33 * pkin(5) - t30 * pkin(13) + t5;
t2 = t67 * t3 + t89 * t4;
t1 = t89 * t3 - t67 * t4;
t8 = [0, 0, 0, 0, 0, t75 / 0.2e1, 0, 0, 0, 0, t71 ^ 2 * t78, t71 * t82, t61 * t80, t74 ^ 2 * t78, t61 * t79, t61 ^ 2 / 0.2e1, pkin(1) * t82 + t53 * t61, -pkin(1) * t71 * t87 - t54 * t61 (-t53 * t71 + t54 * t74) * t84, t54 ^ 2 / 0.2e1 + t53 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t78, t46 ^ 2 / 0.2e1, -t46 * t44, -t46 * t51, t44 ^ 2 / 0.2e1, t44 * t51, t51 ^ 2 / 0.2e1, -t25 * t51 + t31 * t44, t26 * t51 + t31 * t46, -t25 * t46 - t26 * t44, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t36 ^ 2 / 0.2e1, -t34 * t36, t42 * t36, t34 ^ 2 / 0.2e1, -t34 * t42, t42 ^ 2 / 0.2e1, t11 * t42 + t23 * t34, -t12 * t42 + t23 * t36, -t11 * t36 - t12 * t34, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t33 * t30, t28 ^ 2 / 0.2e1, -t33 * t28, t33 ^ 2 / 0.2e1, t9 * t28 + t5 * t33, t9 * t30 - t6 * t33, -t6 * t28 - t5 * t30, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t32, t16 ^ 2 / 0.2e1, -t16 * t32, t32 ^ 2 / 0.2e1, t1 * t32 + t7 * t16, t7 * t18 - t2 * t32, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
