% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:46:28
% EndTime: 2019-03-09 19:46:28
% DurationCPUTime: 0.26s
% Computational Cost: add. (1770->72), mult. (4079->167), div. (0->0), fcn. (3255->12), ass. (0->56)
t65 = cos(qJ(2));
t63 = sin(qJ(2));
t58 = sin(pkin(6));
t74 = qJD(1) * t58;
t69 = t63 * t74;
t73 = cos(pkin(6)) * qJD(1);
t70 = pkin(1) * t73;
t46 = -pkin(8) * t69 + t65 * t70;
t55 = qJD(2) + t73;
t38 = -t55 * pkin(2) - t46;
t62 = sin(qJ(3));
t64 = cos(qJ(3));
t43 = -t64 * t55 + t62 * t69;
t45 = t62 * t55 + t64 * t69;
t25 = t43 * pkin(3) - t45 * qJ(4) + t38;
t68 = t65 * t74;
t47 = pkin(8) * t68 + t63 * t70;
t39 = t55 * pkin(9) + t47;
t41 = (-pkin(2) * t65 - pkin(9) * t63 - pkin(1)) * t74;
t30 = t64 * t39 + t62 * t41;
t50 = -qJD(3) + t68;
t28 = -t50 * qJ(4) + t30;
t57 = sin(pkin(12));
t75 = cos(pkin(12));
t17 = t57 * t25 + t75 * t28;
t32 = t57 * t45 + t75 * t50;
t12 = -t32 * pkin(10) + t17;
t61 = sin(qJ(5));
t78 = cos(qJ(5));
t16 = t75 * t25 - t57 * t28;
t34 = t75 * t45 - t57 * t50;
t9 = t43 * pkin(4) - t34 * pkin(10) + t16;
t6 = t78 * t12 + t61 * t9;
t77 = cos(qJ(6));
t66 = qJD(1) ^ 2;
t76 = t58 ^ 2 * t66;
t72 = t43 ^ 2 / 0.2e1;
t71 = t65 * t76;
t67 = t76 / 0.2e1;
t5 = -t61 * t12 + t78 * t9;
t29 = -t62 * t39 + t64 * t41;
t42 = qJD(5) + t43;
t27 = t50 * pkin(3) + qJD(4) - t29;
t18 = t32 * pkin(4) + t27;
t60 = sin(qJ(6));
t40 = qJD(6) + t42;
t22 = -t61 * t32 + t78 * t34;
t20 = t78 * t32 + t61 * t34;
t15 = -t60 * t20 + t77 * t22;
t13 = t77 * t20 + t60 * t22;
t10 = t20 * pkin(5) + t18;
t4 = -t20 * pkin(11) + t6;
t3 = t42 * pkin(5) - t22 * pkin(11) + t5;
t2 = t60 * t3 + t77 * t4;
t1 = t77 * t3 - t60 * t4;
t7 = [0, 0, 0, 0, 0, t66 / 0.2e1, 0, 0, 0, 0, t63 ^ 2 * t67, t63 * t71, t55 * t69, t65 ^ 2 * t67, t55 * t68, t55 ^ 2 / 0.2e1, pkin(1) * t71 + t46 * t55, -pkin(1) * t63 * t76 - t47 * t55 (-t46 * t63 + t47 * t65) * t74, t47 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t67, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t50, t72, t43 * t50, t50 ^ 2 / 0.2e1, -t29 * t50 + t38 * t43, t30 * t50 + t38 * t45, -t29 * t45 - t30 * t43, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t43, t32 ^ 2 / 0.2e1, -t32 * t43, t72, t16 * t43 + t27 * t32, -t17 * t43 + t27 * t34, -t16 * t34 - t17 * t32, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t42, t20 ^ 2 / 0.2e1, -t20 * t42, t42 ^ 2 / 0.2e1, t18 * t20 + t5 * t42, t18 * t22 - t6 * t42, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t40, t13 ^ 2 / 0.2e1, -t13 * t40, t40 ^ 2 / 0.2e1, t1 * t40 + t10 * t13, t10 * t15 - t2 * t40, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t7;
