% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:59:59
% EndTime: 2019-03-09 00:59:59
% DurationCPUTime: 0.26s
% Computational Cost: add. (988->64), mult. (2373->160), div. (0->0), fcn. (1925->14), ass. (0->57)
t61 = cos(qJ(2));
t52 = sin(pkin(6));
t72 = qJD(1) * t52;
t40 = qJD(2) * pkin(2) + t61 * t72;
t51 = sin(pkin(7));
t53 = cos(pkin(7));
t54 = cos(pkin(6));
t71 = qJD(1) * t54;
t78 = t40 * t53 + t51 * t71;
t59 = sin(qJ(2));
t70 = qJD(2) * t51;
t39 = pkin(9) * t70 + t59 * t72;
t58 = sin(qJ(3));
t60 = cos(qJ(3));
t25 = -t58 * t39 + t78 * t60;
t26 = t60 * t39 + t78 * t58;
t48 = t53 * qJD(2) + qJD(3);
t20 = t48 * pkin(10) + t26;
t47 = t53 * t71;
t29 = t47 + (-t40 + (-pkin(3) * t60 - pkin(10) * t58) * qJD(2)) * t51;
t57 = sin(qJ(4));
t77 = cos(qJ(4));
t13 = t77 * t20 + t57 * t29;
t68 = t58 * t70;
t32 = -t77 * t48 + t57 * t68;
t11 = -t32 * pkin(11) + t13;
t56 = sin(qJ(5));
t76 = cos(qJ(5));
t12 = -t57 * t20 + t77 * t29;
t34 = t57 * t48 + t77 * t68;
t67 = t60 * t70;
t45 = -qJD(4) + t67;
t9 = -t45 * pkin(4) - t34 * pkin(11) + t12;
t6 = t76 * t11 + t56 * t9;
t75 = cos(qJ(6));
t62 = qJD(2) ^ 2;
t73 = t51 ^ 2 * t62;
t66 = t73 / 0.2e1;
t65 = qJD(2) * t72;
t22 = t76 * t32 + t56 * t34;
t5 = -t56 * t11 + t76 * t9;
t19 = -t48 * pkin(3) - t25;
t14 = t32 * pkin(4) + t19;
t63 = qJD(1) ^ 2;
t55 = sin(qJ(6));
t41 = -qJD(5) + t45;
t31 = -t51 * t40 + t47;
t24 = -t56 * t32 + t76 * t34;
t21 = qJD(6) + t22;
t17 = t75 * t24 - t55 * t41;
t15 = t55 * t24 + t75 * t41;
t7 = t22 * pkin(5) - t24 * pkin(12) + t14;
t4 = -t41 * pkin(12) + t6;
t3 = t41 * pkin(5) - t5;
t2 = t75 * t4 + t55 * t7;
t1 = -t55 * t4 + t75 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t63 / 0.2e1, 0, 0, 0, 0, 0, t62 / 0.2e1, t61 * t65, -t59 * t65, 0 (t54 ^ 2 / 0.2e1 + (t59 ^ 2 / 0.2e1 + t61 ^ 2 / 0.2e1) * t52 ^ 2) * t63, t58 ^ 2 * t66, t58 * t60 * t73, t48 * t68, t60 ^ 2 * t66, t48 * t67, t48 ^ 2 / 0.2e1, t25 * t48 - t31 * t67, -t26 * t48 + t31 * t68 (-t25 * t58 + t26 * t60) * t70, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, -t34 * t45, t32 ^ 2 / 0.2e1, t32 * t45, t45 ^ 2 / 0.2e1, -t12 * t45 + t19 * t32, t13 * t45 + t19 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, -t24 * t41, t22 ^ 2 / 0.2e1, t22 * t41, t41 ^ 2 / 0.2e1, t14 * t22 - t5 * t41, t14 * t24 + t6 * t41, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t21, t15 ^ 2 / 0.2e1, -t15 * t21, t21 ^ 2 / 0.2e1, t1 * t21 + t3 * t15, t3 * t17 - t2 * t21, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
