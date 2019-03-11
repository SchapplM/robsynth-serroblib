% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRP10
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:36:54
% EndTime: 2019-03-09 17:36:54
% DurationCPUTime: 0.23s
% Computational Cost: add. (1272->66), mult. (2974->147), div. (0->0), fcn. (2321->10), ass. (0->55)
t55 = sin(qJ(5));
t58 = cos(qJ(2));
t57 = sin(qJ(2));
t53 = sin(pkin(6));
t68 = qJD(1) * t53;
t62 = t57 * t68;
t67 = cos(pkin(6)) * qJD(1);
t63 = pkin(1) * t67;
t41 = -pkin(8) * t62 + t58 * t63;
t50 = qJD(2) + t67;
t33 = -t50 * pkin(2) - t41;
t56 = sin(qJ(3));
t74 = cos(qJ(3));
t38 = -t74 * t50 + t56 * t62;
t40 = t56 * t50 + t74 * t62;
t20 = t38 * pkin(3) - t40 * qJ(4) + t33;
t61 = t58 * t68;
t42 = pkin(8) * t61 + t57 * t63;
t34 = t50 * pkin(9) + t42;
t36 = (-pkin(2) * t58 - pkin(9) * t57 - pkin(1)) * t68;
t25 = t74 * t34 + t56 * t36;
t45 = -qJD(3) + t61;
t23 = -t45 * qJ(4) + t25;
t52 = sin(pkin(11));
t69 = cos(pkin(11));
t10 = t69 * t20 - t52 * t23;
t29 = t69 * t40 - t52 * t45;
t7 = t38 * pkin(4) - t29 * pkin(10) + t10;
t73 = cos(qJ(5));
t11 = t52 * t20 + t69 * t23;
t27 = t52 * t40 + t69 * t45;
t9 = -t27 * pkin(10) + t11;
t4 = t55 * t7 + t73 * t9;
t15 = t73 * t27 + t55 * t29;
t17 = -t55 * t27 + t73 * t29;
t72 = t17 * t15;
t37 = qJD(5) + t38;
t71 = t37 * t15;
t59 = qJD(1) ^ 2;
t70 = t53 ^ 2 * t59;
t66 = t15 ^ 2 / 0.2e1;
t65 = t38 ^ 2 / 0.2e1;
t64 = t58 * t70;
t60 = t70 / 0.2e1;
t24 = -t56 * t34 + t74 * t36;
t3 = -t55 * t9 + t73 * t7;
t22 = t45 * pkin(3) + qJD(4) - t24;
t12 = t27 * pkin(4) + t22;
t35 = t37 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t13 = t17 * t37;
t5 = t15 * pkin(5) - t17 * qJ(6) + t12;
t2 = t37 * qJ(6) + t4;
t1 = -t37 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t59 / 0.2e1, 0, 0, 0, 0, t57 ^ 2 * t60, t57 * t64, t50 * t62, t58 ^ 2 * t60, t50 * t61, t50 ^ 2 / 0.2e1, pkin(1) * t64 + t41 * t50, -pkin(1) * t57 * t70 - t42 * t50 (-t41 * t57 + t42 * t58) * t68, t42 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t60, t40 ^ 2 / 0.2e1, -t40 * t38, -t40 * t45, t65, t38 * t45, t45 ^ 2 / 0.2e1, -t24 * t45 + t33 * t38, t25 * t45 + t33 * t40, -t24 * t40 - t25 * t38, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t38, t27 ^ 2 / 0.2e1, -t27 * t38, t65, t10 * t38 + t22 * t27, -t11 * t38 + t22 * t29, -t10 * t29 - t11 * t27, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t14, -t72, t13, t66, -t71, t35, t12 * t15 + t3 * t37, t12 * t17 - t4 * t37, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t14, t13, t72, t35, t71, t66, -t1 * t37 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t37, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
