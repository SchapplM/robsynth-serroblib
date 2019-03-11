% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:35
% EndTime: 2019-03-09 12:48:35
% DurationCPUTime: 0.16s
% Computational Cost: add. (539->62), mult. (1139->128), div. (0->0), fcn. (670->6), ass. (0->52)
t53 = qJD(1) ^ 2;
t65 = t53 / 0.2e1;
t64 = -pkin(2) - pkin(8);
t52 = cos(qJ(2));
t50 = sin(qJ(2));
t55 = -qJ(3) * t50 - pkin(1);
t25 = (t64 * t52 + t55) * qJD(1);
t60 = t50 * qJD(1);
t59 = pkin(7) * t60 + qJD(3);
t26 = pkin(3) * t60 + t64 * qJD(2) + t59;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t13 = t51 * t25 + t49 * t26;
t61 = qJD(1) * t52;
t30 = t49 * qJD(2) + t51 * t61;
t11 = -t30 * pkin(9) + t13;
t48 = sin(qJ(5));
t63 = cos(qJ(5));
t12 = -t49 * t25 + t51 * t26;
t32 = t51 * qJD(2) - t49 * t61;
t38 = qJD(4) + t60;
t7 = t38 * pkin(4) - t32 * pkin(9) + t12;
t4 = t63 * t11 + t48 * t7;
t62 = t52 * t53;
t35 = -pkin(7) * t61 - qJD(2) * qJ(3);
t58 = qJD(1) * qJD(2);
t28 = pkin(3) * t61 - t35;
t3 = -t48 * t11 + t63 * t7;
t57 = t50 * t58;
t56 = t52 * t58;
t21 = t30 * pkin(4) + t28;
t47 = t52 ^ 2;
t46 = t50 ^ 2;
t44 = qJD(2) ^ 2 / 0.2e1;
t40 = t47 * t65;
t39 = t46 * t65;
t37 = t50 * t62;
t36 = qJD(5) + t38;
t34 = t36 ^ 2 / 0.2e1;
t33 = -qJD(2) * pkin(2) + t59;
t29 = (-pkin(2) * t52 + t55) * qJD(1);
t20 = -t48 * t30 + t63 * t32;
t18 = t63 * t30 + t48 * t32;
t17 = t20 ^ 2 / 0.2e1;
t16 = t18 ^ 2 / 0.2e1;
t15 = t36 * t20;
t14 = t36 * t18;
t9 = t18 * pkin(5) + qJD(6) + t21;
t8 = t20 * t18;
t2 = -t18 * qJ(6) + t4;
t1 = t36 * pkin(5) - t20 * qJ(6) + t3;
t5 = [0, 0, 0, 0, 0, t65, 0, 0, 0, 0, t39, t37, t57, t40, t56, t44, pkin(1) * t62 - pkin(7) * t57, -t53 * pkin(1) * t50 - pkin(7) * t56 (t46 + t47) * t53 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t47 / 0.2e1 + t46 / 0.2e1) * pkin(7) ^ 2) * t53, t44, -t57, -t56, t39, t37, t40 (t33 * t50 - t35 * t52) * qJD(1), t33 * qJD(2) + t29 * t61, -t35 * qJD(2) - t29 * t60, t29 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t30, t32 * t38, t30 ^ 2 / 0.2e1, -t30 * t38, t38 ^ 2 / 0.2e1, t12 * t38 + t28 * t30, -t13 * t38 + t28 * t32, -t12 * t32 - t13 * t30, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t17, -t8, t15, t16, -t14, t34, t21 * t18 + t3 * t36, t21 * t20 - t4 * t36, -t4 * t18 - t3 * t20, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t17, -t8, t15, t16, -t14, t34, t1 * t36 + t9 * t18, -t2 * t36 + t9 * t20, -t1 * t20 - t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg  = t5;
