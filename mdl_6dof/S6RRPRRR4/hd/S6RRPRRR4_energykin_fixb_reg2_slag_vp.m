% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:34:48
% EndTime: 2019-03-09 13:34:48
% DurationCPUTime: 0.26s
% Computational Cost: add. (1412->73), mult. (3882->168), div. (0->0), fcn. (3100->12), ass. (0->57)
t64 = cos(qJ(2));
t71 = cos(pkin(6)) * qJD(1);
t69 = pkin(1) * t71;
t53 = t64 * t69;
t54 = qJD(2) + t71;
t63 = sin(qJ(2));
t57 = sin(pkin(6));
t72 = qJD(1) * t57;
t68 = t63 * t72;
t35 = t54 * pkin(2) + t53 + (-pkin(8) - qJ(3)) * t68;
t67 = t64 * t72;
t46 = pkin(8) * t67 + t63 * t69;
t39 = qJ(3) * t67 + t46;
t56 = sin(pkin(12));
t58 = cos(pkin(12));
t26 = t56 * t35 + t58 * t39;
t24 = t54 * pkin(9) + t26;
t42 = t56 * t68 - t58 * t67;
t44 = (t56 * t64 + t58 * t63) * t72;
t47 = qJD(3) + (-pkin(2) * t64 - pkin(1)) * t72;
t29 = t42 * pkin(3) - t44 * pkin(9) + t47;
t62 = sin(qJ(4));
t76 = cos(qJ(4));
t13 = t76 * t24 + t62 * t29;
t32 = t62 * t44 - t76 * t54;
t11 = -t32 * pkin(10) + t13;
t61 = sin(qJ(5));
t75 = cos(qJ(5));
t12 = -t62 * t24 + t76 * t29;
t34 = t76 * t44 + t62 * t54;
t41 = qJD(4) + t42;
t9 = t41 * pkin(4) - t34 * pkin(10) + t12;
t6 = t75 * t11 + t61 * t9;
t74 = cos(qJ(6));
t65 = qJD(1) ^ 2;
t73 = t57 ^ 2 * t65;
t70 = t64 * t73;
t66 = t73 / 0.2e1;
t20 = t75 * t32 + t61 * t34;
t25 = t58 * t35 - t56 * t39;
t23 = -t54 * pkin(3) - t25;
t5 = -t61 * t11 + t75 * t9;
t14 = t32 * pkin(4) + t23;
t60 = sin(qJ(6));
t50 = t54 ^ 2 / 0.2e1;
t45 = -pkin(8) * t68 + t53;
t40 = qJD(5) + t41;
t22 = -t61 * t32 + t75 * t34;
t19 = qJD(6) + t20;
t17 = t74 * t22 + t60 * t40;
t15 = t60 * t22 - t74 * t40;
t7 = t20 * pkin(5) - t22 * pkin(11) + t14;
t4 = t40 * pkin(11) + t6;
t3 = -t40 * pkin(5) - t5;
t2 = t74 * t4 + t60 * t7;
t1 = -t60 * t4 + t74 * t7;
t8 = [0, 0, 0, 0, 0, t65 / 0.2e1, 0, 0, 0, 0, t63 ^ 2 * t66, t63 * t70, t54 * t68, t64 ^ 2 * t66, t54 * t67, t50, pkin(1) * t70 + t45 * t54, -pkin(1) * t63 * t73 - t46 * t54 (-t45 * t63 + t46 * t64) * t72, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t66, t44 ^ 2 / 0.2e1, -t44 * t42, t44 * t54, t42 ^ 2 / 0.2e1, -t42 * t54, t50, t25 * t54 + t47 * t42, -t26 * t54 + t47 * t44, -t25 * t44 - t26 * t42, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t41, t32 ^ 2 / 0.2e1, -t32 * t41, t41 ^ 2 / 0.2e1, t12 * t41 + t23 * t32, -t13 * t41 + t23 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t40, t20 ^ 2 / 0.2e1, -t20 * t40, t40 ^ 2 / 0.2e1, t14 * t20 + t5 * t40, t14 * t22 - t6 * t40, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t19, t15 ^ 2 / 0.2e1, -t15 * t19, t19 ^ 2 / 0.2e1, t1 * t19 + t3 * t15, t3 * t17 - t2 * t19, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
