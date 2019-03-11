% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR14_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:36:55
% EndTime: 2019-03-09 11:36:55
% DurationCPUTime: 0.21s
% Computational Cost: add. (610->68), mult. (1467->136), div. (0->0), fcn. (1013->8), ass. (0->61)
t73 = -pkin(2) - pkin(9);
t72 = pkin(4) + pkin(10);
t71 = cos(qJ(6));
t70 = sin(qJ(4));
t46 = cos(pkin(6));
t64 = t46 * qJD(1);
t43 = qJD(2) + t64;
t49 = cos(qJ(4));
t50 = cos(qJ(2));
t45 = sin(pkin(6));
t65 = qJD(1) * t45;
t59 = t50 * t65;
t25 = t70 * t43 + t49 * t59;
t27 = t49 * t43 - t70 * t59;
t69 = t27 * t25;
t48 = sin(qJ(2));
t58 = t48 * t65;
t34 = qJD(4) + t58;
t68 = t27 * t34;
t67 = t34 * t25;
t51 = qJD(1) ^ 2;
t66 = t45 ^ 2 * t51;
t39 = pkin(8) * t58;
t12 = qJD(3) + t39 + t73 * t43 + (-pkin(1) * t46 * t50 + pkin(3) * t45 * t48) * qJD(1);
t56 = -qJ(3) * t48 - pkin(1);
t19 = (t73 * t50 + t56) * t65;
t9 = t70 * t12 + t49 * t19;
t60 = pkin(1) * t64;
t30 = pkin(8) * t59 + t48 * t60;
t63 = t25 ^ 2 / 0.2e1;
t62 = t27 ^ 2 / 0.2e1;
t61 = t50 * t66;
t21 = -t43 * qJ(3) - t30;
t57 = t66 / 0.2e1;
t8 = t49 * t12 - t70 * t19;
t18 = pkin(3) * t59 - t21;
t55 = t43 * t58;
t54 = t43 * t59;
t7 = -t34 * qJ(5) - t9;
t53 = qJD(5) - t8;
t29 = t50 * t60 - t39;
t52 = -t27 * qJ(5) + t18;
t47 = sin(qJ(6));
t37 = t50 ^ 2 * t57;
t36 = t48 ^ 2 * t57;
t35 = t43 ^ 2 / 0.2e1;
t33 = t48 * t61;
t31 = t34 ^ 2 / 0.2e1;
t24 = qJD(6) + t27;
t23 = (-pkin(2) * t50 + t56) * t65;
t20 = -t43 * pkin(2) + qJD(3) - t29;
t15 = t47 * t25 + t71 * t34;
t13 = -t71 * t25 + t47 * t34;
t10 = t25 * pkin(4) + t52;
t6 = -t34 * pkin(4) + t53;
t5 = t72 * t25 + t52;
t4 = -t25 * pkin(5) - t7;
t3 = t27 * pkin(5) - t72 * t34 + t53;
t2 = t47 * t3 + t71 * t5;
t1 = t71 * t3 - t47 * t5;
t11 = [0, 0, 0, 0, 0, t51 / 0.2e1, 0, 0, 0, 0, t36, t33, t55, t37, t54, t35, pkin(1) * t61 + t29 * t43, -pkin(1) * t48 * t66 - t30 * t43 (-t29 * t48 + t30 * t50) * t65, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t57, t35, -t55, -t54, t36, t33, t37 (t20 * t48 - t21 * t50) * t65, t20 * t43 + t23 * t59, -t21 * t43 - t23 * t58, t23 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t62, -t69, t68, t63, -t67, t31, t18 * t25 + t8 * t34, t18 * t27 - t9 * t34, -t9 * t25 - t8 * t27, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t31, -t68, t67, t62, -t69, t63, t7 * t25 + t6 * t27, -t10 * t25 + t6 * t34, -t10 * t27 - t7 * t34, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t24, t13 ^ 2 / 0.2e1, -t13 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t4 * t13, t4 * t15 - t2 * t24, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
