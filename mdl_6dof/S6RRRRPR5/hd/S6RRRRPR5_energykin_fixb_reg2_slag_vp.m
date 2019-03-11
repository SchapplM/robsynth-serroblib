% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:16:15
% EndTime: 2019-03-09 22:16:16
% DurationCPUTime: 0.22s
% Computational Cost: add. (769->63), mult. (1651->143), div. (0->0), fcn. (1145->8), ass. (0->56)
t53 = qJD(1) ^ 2;
t72 = t53 / 0.2e1;
t71 = -pkin(4) - pkin(5);
t70 = -pkin(8) - pkin(7);
t69 = cos(qJ(3));
t68 = cos(qJ(4));
t67 = cos(qJ(6));
t50 = sin(qJ(3));
t51 = sin(qJ(2));
t52 = cos(qJ(2));
t37 = (t50 * t52 + t69 * t51) * qJD(1);
t45 = qJD(2) + qJD(3);
t49 = sin(qJ(4));
t25 = t49 * t37 - t68 * t45;
t27 = t68 * t37 + t49 * t45;
t66 = t27 * t25;
t62 = qJD(1) * t52;
t63 = qJD(1) * t51;
t35 = t50 * t63 - t69 * t62;
t34 = qJD(4) + t35;
t65 = t34 * t25;
t64 = t52 * t53;
t42 = (-pkin(2) * t52 - pkin(1)) * qJD(1);
t16 = t35 * pkin(3) - t37 * pkin(9) + t42;
t40 = qJD(2) * pkin(2) + t70 * t63;
t41 = t70 * t62;
t22 = t50 * t40 - t69 * t41;
t20 = t45 * pkin(9) + t22;
t10 = t49 * t16 + t68 * t20;
t21 = t69 * t40 + t50 * t41;
t61 = t25 ^ 2 / 0.2e1;
t60 = qJD(1) * qJD(2);
t7 = t34 * qJ(5) + t10;
t59 = t51 * t60;
t58 = t52 * t60;
t57 = t45 * pkin(3) + t21;
t9 = t68 * t16 - t49 * t20;
t56 = qJD(5) - t9;
t55 = t27 * qJ(5) + t57;
t48 = sin(qJ(6));
t47 = t52 ^ 2;
t46 = t51 ^ 2;
t30 = -qJD(6) + t34;
t28 = t34 ^ 2 / 0.2e1;
t24 = t27 ^ 2 / 0.2e1;
t15 = t27 * t34;
t13 = t48 * t25 + t67 * t27;
t11 = -t67 * t25 + t48 * t27;
t8 = t25 * pkin(4) - t55;
t6 = -t34 * pkin(4) + t56;
t5 = t71 * t25 + t55;
t4 = t25 * pkin(10) + t7;
t3 = -t27 * pkin(10) + t71 * t34 + t56;
t2 = t48 * t3 + t67 * t4;
t1 = t67 * t3 - t48 * t4;
t12 = [0, 0, 0, 0, 0, t72, 0, 0, 0, 0, t46 * t72, t51 * t64, t59, t47 * t72, t58, qJD(2) ^ 2 / 0.2e1, pkin(1) * t64 - pkin(7) * t59, -t53 * pkin(1) * t51 - pkin(7) * t58 (t46 + t47) * t53 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t47 / 0.2e1 + t46 / 0.2e1) * pkin(7) ^ 2) * t53, t37 ^ 2 / 0.2e1, -t37 * t35, t37 * t45, t35 ^ 2 / 0.2e1, -t35 * t45, t45 ^ 2 / 0.2e1, t21 * t45 + t42 * t35, -t22 * t45 + t42 * t37, -t21 * t37 - t22 * t35, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1, t24, -t66, t15, t61, -t65, t28, -t25 * t57 + t9 * t34, -t10 * t34 - t27 * t57, -t10 * t25 - t9 * t27, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t57 ^ 2 / 0.2e1, t24, t15, t66, t28, t65, t61, t8 * t25 - t6 * t34, -t7 * t25 + t6 * t27, -t8 * t27 + t7 * t34, t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t11 * t13, -t30 * t13, t11 ^ 2 / 0.2e1, t11 * t30, t30 ^ 2 / 0.2e1, -t1 * t30 + t5 * t11, t5 * t13 + t2 * t30, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
