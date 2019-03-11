% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:51:24
% EndTime: 2019-03-10 03:51:24
% DurationCPUTime: 0.24s
% Computational Cost: add. (1409->69), mult. (3064->164), div. (0->0), fcn. (2249->10), ass. (0->53)
t63 = qJD(1) ^ 2;
t75 = t63 / 0.2e1;
t61 = sin(qJ(2));
t62 = cos(qJ(2));
t38 = (-pkin(2) * t62 - pkin(8) * t61 - pkin(1)) * qJD(1);
t68 = t62 * qJD(1);
t48 = pkin(7) * t68 + qJD(2) * pkin(8);
t60 = sin(qJ(3));
t74 = cos(qJ(3));
t32 = t74 * t38 - t60 * t48;
t69 = qJD(1) * t61;
t41 = t60 * qJD(2) + t74 * t69;
t51 = -qJD(3) + t68;
t25 = -t51 * pkin(3) - t41 * pkin(9) + t32;
t33 = t60 * t38 + t74 * t48;
t39 = -t74 * qJD(2) + t60 * t69;
t27 = -t39 * pkin(9) + t33;
t59 = sin(qJ(4));
t73 = cos(qJ(4));
t16 = t73 * t25 - t59 * t27;
t31 = -t59 * t39 + t73 * t41;
t49 = -qJD(4) + t51;
t12 = -t49 * pkin(4) - t31 * pkin(10) + t16;
t17 = t59 * t25 + t73 * t27;
t29 = t73 * t39 + t59 * t41;
t14 = -t29 * pkin(10) + t17;
t58 = sin(qJ(5));
t72 = cos(qJ(5));
t6 = t58 * t12 + t72 * t14;
t71 = cos(qJ(6));
t70 = t62 * t63;
t67 = qJD(1) * qJD(2);
t5 = t72 * t12 - t58 * t14;
t66 = t61 * t67;
t65 = t62 * t67;
t47 = -qJD(2) * pkin(2) + pkin(7) * t69;
t34 = t39 * pkin(3) + t47;
t45 = -qJD(5) + t49;
t22 = t29 * pkin(4) + t34;
t57 = sin(qJ(6));
t56 = t62 ^ 2;
t55 = t61 ^ 2;
t43 = -qJD(6) + t45;
t21 = -t58 * t29 + t72 * t31;
t19 = t72 * t29 + t58 * t31;
t15 = t19 * pkin(5) + t22;
t11 = -t57 * t19 + t71 * t21;
t9 = t71 * t19 + t57 * t21;
t4 = -t19 * pkin(11) + t6;
t3 = -t45 * pkin(5) - t21 * pkin(11) + t5;
t2 = t57 * t3 + t71 * t4;
t1 = t71 * t3 - t57 * t4;
t7 = [0, 0, 0, 0, 0, t75, 0, 0, 0, 0, t55 * t75, t61 * t70, t66, t56 * t75, t65, qJD(2) ^ 2 / 0.2e1, pkin(1) * t70 - pkin(7) * t66, -t63 * pkin(1) * t61 - pkin(7) * t65 (t55 + t56) * t63 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t56 / 0.2e1 + t55 / 0.2e1) * pkin(7) ^ 2) * t63, t41 ^ 2 / 0.2e1, -t41 * t39, -t41 * t51, t39 ^ 2 / 0.2e1, t39 * t51, t51 ^ 2 / 0.2e1, -t32 * t51 + t47 * t39, t33 * t51 + t47 * t41, -t32 * t41 - t33 * t39, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t49, t29 ^ 2 / 0.2e1, t29 * t49, t49 ^ 2 / 0.2e1, -t16 * t49 + t34 * t29, t17 * t49 + t34 * t31, -t16 * t31 - t17 * t29, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, -t21 * t45, t19 ^ 2 / 0.2e1, t19 * t45, t45 ^ 2 / 0.2e1, t22 * t19 - t5 * t45, t22 * t21 + t6 * t45, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t9, -t11 * t43, t9 ^ 2 / 0.2e1, t9 * t43, t43 ^ 2 / 0.2e1, -t1 * t43 + t15 * t9, t15 * t11 + t2 * t43, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1;];
T_reg  = t7;
