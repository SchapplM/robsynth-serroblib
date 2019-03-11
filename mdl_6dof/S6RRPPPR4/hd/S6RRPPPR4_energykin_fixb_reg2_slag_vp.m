% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:43
% EndTime: 2019-03-09 08:19:43
% DurationCPUTime: 0.16s
% Computational Cost: add. (409->61), mult. (914->127), div. (0->0), fcn. (494->6), ass. (0->53)
t48 = qJD(1) ^ 2;
t66 = t48 / 0.2e1;
t65 = -pkin(4) - pkin(5);
t64 = cos(qJ(6));
t43 = sin(pkin(9));
t44 = cos(pkin(9));
t47 = cos(qJ(2));
t60 = qJD(1) * t47;
t24 = t43 * qJD(2) + t44 * t60;
t26 = t44 * qJD(2) - t43 * t60;
t63 = t26 * t24;
t62 = t47 * t48;
t61 = -pkin(2) - qJ(4);
t46 = sin(qJ(2));
t52 = -qJ(3) * t46 - pkin(1);
t17 = (t61 * t47 + t52) * qJD(1);
t59 = t46 * qJD(1);
t57 = pkin(7) * t59 + qJD(3);
t18 = pkin(3) * t59 + t61 * qJD(2) + t57;
t9 = t44 * t17 + t43 * t18;
t28 = -pkin(7) * t60 - qJD(2) * qJ(3);
t58 = t24 ^ 2 / 0.2e1;
t56 = qJD(1) * qJD(2);
t7 = qJ(5) * t59 + t9;
t55 = t24 * t59;
t54 = t46 * t56;
t53 = t47 * t56;
t8 = -t43 * t17 + t44 * t18;
t21 = pkin(3) * t60 + qJD(4) - t28;
t51 = qJD(5) - t8;
t50 = t26 * qJ(5) - t21;
t45 = sin(qJ(6));
t42 = t47 ^ 2;
t41 = t46 ^ 2;
t39 = qJD(2) ^ 2 / 0.2e1;
t34 = t42 * t66;
t33 = t41 * t66;
t31 = -qJD(6) + t59;
t30 = t46 * t62;
t27 = -qJD(2) * pkin(2) + t57;
t23 = (-pkin(2) * t47 + t52) * qJD(1);
t22 = t26 ^ 2 / 0.2e1;
t19 = t26 * t59;
t13 = t45 * t24 + t64 * t26;
t11 = -t64 * t24 + t45 * t26;
t10 = t24 * pkin(4) - t50;
t6 = -pkin(4) * t59 + t51;
t5 = t65 * t24 + t50;
t4 = t24 * pkin(8) + t7;
t3 = -t26 * pkin(8) + t65 * t59 + t51;
t2 = t45 * t3 + t64 * t4;
t1 = t64 * t3 - t45 * t4;
t12 = [0, 0, 0, 0, 0, t66, 0, 0, 0, 0, t33, t30, t54, t34, t53, t39, pkin(1) * t62 - pkin(7) * t54, -t48 * pkin(1) * t46 - pkin(7) * t53 (t41 + t42) * t48 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * pkin(7) ^ 2) * t48, t39, -t54, -t53, t33, t30, t34 (t27 * t46 - t28 * t47) * qJD(1), t27 * qJD(2) + t23 * t60, -t28 * qJD(2) - t23 * t59, t23 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t22, -t63, t19, t58, -t55, t33, t21 * t24 + t8 * t59, t21 * t26 - t9 * t59, -t9 * t24 - t8 * t26, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t22, t19, t63, t33, t55, t58, t10 * t24 - t6 * t59, -t7 * t24 + t6 * t26, -t10 * t26 + t7 * t59, t7 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, -t13 * t31, t11 ^ 2 / 0.2e1, t11 * t31, t31 ^ 2 / 0.2e1, -t1 * t31 + t5 * t11, t5 * t13 + t2 * t31, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
