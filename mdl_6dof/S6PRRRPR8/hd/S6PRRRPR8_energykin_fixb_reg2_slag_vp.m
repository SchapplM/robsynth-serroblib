% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:53:19
% EndTime: 2019-03-08 23:53:19
% DurationCPUTime: 0.23s
% Computational Cost: add. (622->61), mult. (1532->140), div. (0->0), fcn. (1195->12), ass. (0->59)
t53 = cos(qJ(2));
t44 = sin(pkin(6));
t68 = qJD(1) * t44;
t34 = qJD(2) * pkin(2) + t53 * t68;
t43 = sin(pkin(7));
t45 = cos(pkin(7));
t46 = cos(pkin(6));
t67 = qJD(1) * t46;
t76 = t34 * t45 + t43 * t67;
t50 = sin(qJ(2));
t66 = qJD(2) * t43;
t32 = pkin(9) * t66 + t50 * t68;
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t15 = -t49 * t32 + t76 * t52;
t75 = pkin(4) + pkin(11);
t74 = cos(qJ(6));
t40 = t45 * qJD(2) + qJD(3);
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t62 = t49 * t66;
t25 = -t51 * t40 + t48 * t62;
t27 = t48 * t40 + t51 * t62;
t73 = t27 * t25;
t61 = t52 * t66;
t37 = -qJD(4) + t61;
t72 = t27 * t37;
t70 = t37 * t25;
t54 = qJD(2) ^ 2;
t69 = t43 ^ 2 * t54;
t16 = t52 * t32 + t76 * t49;
t14 = t40 * pkin(10) + t16;
t39 = t45 * t67;
t18 = t39 + (-t34 + (-pkin(3) * t52 - pkin(10) * t49) * qJD(2)) * t43;
t9 = t51 * t14 + t48 * t18;
t65 = t25 ^ 2 / 0.2e1;
t64 = t27 ^ 2 / 0.2e1;
t60 = t69 / 0.2e1;
t59 = qJD(2) * t68;
t8 = -t48 * t14 + t51 * t18;
t7 = t37 * qJ(5) - t9;
t58 = qJD(5) - t8;
t13 = -t40 * pkin(3) - t15;
t56 = -t27 * qJ(5) + t13;
t55 = qJD(1) ^ 2;
t47 = sin(qJ(6));
t33 = t37 ^ 2 / 0.2e1;
t24 = qJD(6) + t27;
t23 = -t43 * t34 + t39;
t21 = t47 * t25 - t74 * t37;
t19 = -t74 * t25 - t47 * t37;
t10 = t25 * pkin(4) + t56;
t6 = t37 * pkin(4) + t58;
t5 = t75 * t25 + t56;
t4 = -t25 * pkin(5) - t7;
t3 = t27 * pkin(5) + t75 * t37 + t58;
t2 = t47 * t3 + t74 * t5;
t1 = t74 * t3 - t47 * t5;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t55 / 0.2e1, 0, 0, 0, 0, 0, t54 / 0.2e1, t53 * t59, -t50 * t59, 0 (t46 ^ 2 / 0.2e1 + (t50 ^ 2 / 0.2e1 + t53 ^ 2 / 0.2e1) * t44 ^ 2) * t55, t49 ^ 2 * t60, t49 * t52 * t69, t40 * t62, t52 ^ 2 * t60, t40 * t61, t40 ^ 2 / 0.2e1, t15 * t40 - t23 * t61, -t16 * t40 + t23 * t62 (-t15 * t49 + t16 * t52) * t66, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t64, -t73, -t72, t65, t70, t33, t13 * t25 - t8 * t37, t13 * t27 + t9 * t37, -t9 * t25 - t8 * t27, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t33, t72, -t70, t64, -t73, t65, t7 * t25 + t6 * t27, -t10 * t25 - t6 * t37, -t10 * t27 + t7 * t37, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * t24, t19 ^ 2 / 0.2e1, -t19 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t4 * t19, -t2 * t24 + t4 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
