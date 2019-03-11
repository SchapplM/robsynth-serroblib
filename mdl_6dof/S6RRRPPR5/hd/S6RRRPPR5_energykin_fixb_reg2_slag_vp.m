% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:44:12
% EndTime: 2019-03-09 15:44:12
% DurationCPUTime: 0.26s
% Computational Cost: add. (1663->72), mult. (3927->164), div. (0->0), fcn. (3136->12), ass. (0->56)
t63 = sin(qJ(2));
t64 = cos(qJ(2));
t59 = sin(pkin(6));
t73 = qJD(1) * t59;
t67 = t64 * t73;
t72 = cos(pkin(6)) * qJD(1);
t69 = pkin(1) * t72;
t46 = pkin(8) * t67 + t63 * t69;
t55 = qJD(2) + t72;
t39 = pkin(9) * t55 + t46;
t41 = (-pkin(2) * t64 - pkin(9) * t63 - pkin(1)) * t73;
t62 = sin(qJ(3));
t78 = cos(qJ(3));
t28 = -t39 * t62 + t41 * t78;
t68 = t63 * t73;
t44 = t55 * t62 + t68 * t78;
t50 = -qJD(3) + t67;
t20 = -pkin(3) * t50 - qJ(4) * t44 + t28;
t29 = t39 * t78 + t41 * t62;
t42 = -t55 * t78 + t62 * t68;
t23 = -qJ(4) * t42 + t29;
t58 = sin(pkin(11));
t75 = cos(pkin(11));
t12 = t20 * t58 + t23 * t75;
t10 = -qJ(5) * t50 + t12;
t45 = -pkin(8) * t68 + t64 * t69;
t38 = -pkin(2) * t55 - t45;
t31 = pkin(3) * t42 + qJD(4) + t38;
t32 = t42 * t75 + t44 * t58;
t34 = -t42 * t58 + t44 * t75;
t15 = pkin(4) * t32 - qJ(5) * t34 + t31;
t57 = sin(pkin(12));
t74 = cos(pkin(12));
t6 = t10 * t74 + t15 * t57;
t77 = cos(qJ(6));
t65 = qJD(1) ^ 2;
t76 = t59 ^ 2 * t65;
t71 = t32 ^ 2 / 0.2e1;
t70 = t64 * t76;
t66 = t76 / 0.2e1;
t5 = -t10 * t57 + t15 * t74;
t11 = t20 * t75 - t23 * t58;
t9 = pkin(4) * t50 + qJD(5) - t11;
t61 = sin(qJ(6));
t48 = t50 ^ 2 / 0.2e1;
t30 = qJD(6) + t32;
t27 = t34 * t74 - t50 * t57;
t25 = t34 * t57 + t50 * t74;
t18 = -t25 * t61 + t27 * t77;
t16 = t25 * t77 + t27 * t61;
t7 = pkin(5) * t25 + t9;
t4 = -pkin(10) * t25 + t6;
t3 = pkin(5) * t32 - pkin(10) * t27 + t5;
t2 = t3 * t61 + t4 * t77;
t1 = t3 * t77 - t4 * t61;
t8 = [0, 0, 0, 0, 0, t65 / 0.2e1, 0, 0, 0, 0, t63 ^ 2 * t66, t63 * t70, t55 * t68, t64 ^ 2 * t66, t55 * t67, t55 ^ 2 / 0.2e1, pkin(1) * t70 + t45 * t55, -pkin(1) * t63 * t76 - t46 * t55 (-t45 * t63 + t46 * t64) * t73, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t66, t44 ^ 2 / 0.2e1, -t44 * t42, -t44 * t50, t42 ^ 2 / 0.2e1, t42 * t50, t48, -t28 * t50 + t38 * t42, t29 * t50 + t38 * t44, -t28 * t44 - t29 * t42, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, -t34 * t50, t71, t32 * t50, t48, -t11 * t50 + t31 * t32, t12 * t50 + t31 * t34, -t11 * t34 - t12 * t32, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t32, t25 ^ 2 / 0.2e1, -t25 * t32, t71, t25 * t9 + t32 * t5, t27 * t9 - t32 * t6, -t25 * t6 - t27 * t5, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t30, t16 ^ 2 / 0.2e1, -t16 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t16 * t7, t18 * t7 - t2 * t30, -t1 * t18 - t16 * t2, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
