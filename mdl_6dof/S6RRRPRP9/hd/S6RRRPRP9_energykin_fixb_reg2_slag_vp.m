% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:27:13
% EndTime: 2019-03-09 17:27:14
% DurationCPUTime: 0.16s
% Computational Cost: add. (543->60), mult. (1141->124), div. (0->0), fcn. (710->6), ass. (0->52)
t50 = qJD(1) ^ 2;
t67 = t50 / 0.2e1;
t46 = sin(qJ(5));
t65 = cos(qJ(5));
t47 = sin(qJ(3));
t48 = sin(qJ(2));
t59 = qJD(1) * t48;
t66 = cos(qJ(3));
t30 = t47 * qJD(2) + t66 * t59;
t49 = cos(qJ(2));
t58 = t49 * qJD(1);
t39 = -qJD(3) + t58;
t24 = (-pkin(2) * t49 - pkin(8) * t48 - pkin(1)) * qJD(1);
t35 = pkin(7) * t58 + qJD(2) * pkin(8);
t19 = t66 * t24 - t47 * t35;
t52 = qJD(4) - t19;
t7 = -t30 * pkin(9) + (pkin(3) + pkin(4)) * t39 + t52;
t20 = t47 * t24 + t66 * t35;
t13 = -t39 * qJ(4) + t20;
t28 = -t66 * qJD(2) + t47 * t59;
t9 = t28 * pkin(9) + t13;
t4 = t46 * t7 + t65 * t9;
t16 = -t65 * t28 + t46 * t30;
t18 = t46 * t28 + t65 * t30;
t64 = t18 * t16;
t63 = t30 * t28;
t38 = qJD(5) + t39;
t62 = t38 * t16;
t61 = t39 * t28;
t60 = t49 * t50;
t34 = -qJD(2) * pkin(2) + pkin(7) * t59;
t57 = t16 ^ 2 / 0.2e1;
t56 = t28 ^ 2 / 0.2e1;
t55 = qJD(1) * qJD(2);
t54 = t48 * t55;
t53 = t49 * t55;
t14 = t28 * pkin(3) - t30 * qJ(4) + t34;
t3 = -t46 * t9 + t65 * t7;
t10 = -t28 * pkin(4) - t14;
t44 = t49 ^ 2;
t43 = t48 ^ 2;
t36 = t39 ^ 2 / 0.2e1;
t33 = t38 ^ 2 / 0.2e1;
t25 = t30 ^ 2 / 0.2e1;
t21 = t30 * t39;
t15 = t18 ^ 2 / 0.2e1;
t12 = t39 * pkin(3) + t52;
t11 = t18 * t38;
t5 = t16 * pkin(5) - t18 * qJ(6) + t10;
t2 = t38 * qJ(6) + t4;
t1 = -t38 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t43 * t67, t48 * t60, t54, t44 * t67, t53, qJD(2) ^ 2 / 0.2e1, pkin(1) * t60 - pkin(7) * t54, -t50 * pkin(1) * t48 - pkin(7) * t53 (t43 + t44) * t50 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t44 / 0.2e1 + t43 / 0.2e1) * pkin(7) ^ 2) * t50, t25, -t63, -t21, t56, t61, t36, -t19 * t39 + t34 * t28, t20 * t39 + t34 * t30, -t19 * t30 - t20 * t28, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t25, -t21, t63, t36, -t61, t56, t12 * t39 + t14 * t28, t12 * t30 - t13 * t28, -t13 * t39 - t14 * t30, t13 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t15, -t64, t11, t57, -t62, t33, t10 * t16 + t3 * t38, t10 * t18 - t4 * t38, -t4 * t16 - t3 * t18, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t15, t11, t64, t33, t62, t57, -t1 * t38 + t5 * t16, t1 * t18 - t2 * t16, -t5 * t18 + t2 * t38, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
