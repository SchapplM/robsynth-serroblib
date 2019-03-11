% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:05
% EndTime: 2019-03-09 08:12:05
% DurationCPUTime: 0.22s
% Computational Cost: add. (646->65), mult. (1578->139), div. (0->0), fcn. (1072->8), ass. (0->54)
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t48 = sin(qJ(2));
t49 = cos(qJ(2));
t32 = (t45 * t49 + t46 * t48) * qJD(1);
t50 = qJD(1) ^ 2;
t69 = t50 / 0.2e1;
t61 = qJD(1) * t49;
t62 = qJD(1) * t48;
t30 = t45 * t62 - t46 * t61;
t38 = qJD(3) + (-pkin(2) * t49 - pkin(1)) * qJD(1);
t52 = -t32 * qJ(4) + t38;
t65 = pkin(3) + qJ(5);
t11 = t65 * t30 + t52;
t64 = pkin(7) + qJ(3);
t36 = qJD(2) * pkin(2) - t64 * t62;
t37 = t64 * t61;
t20 = t46 * t36 - t45 * t37;
t54 = qJD(4) - t20;
t12 = t32 * pkin(4) - t65 * qJD(2) + t54;
t44 = sin(pkin(10));
t63 = cos(pkin(10));
t6 = t63 * t11 + t44 * t12;
t68 = cos(qJ(6));
t67 = t32 * t30;
t66 = t49 * t50;
t21 = t45 * t36 + t46 * t37;
t60 = qJD(2) * t30;
t59 = t32 * qJD(2);
t58 = t30 ^ 2 / 0.2e1;
t28 = t32 ^ 2 / 0.2e1;
t57 = qJD(1) * qJD(2);
t19 = -qJD(2) * qJ(4) - t21;
t56 = t48 * t57;
t55 = t49 * t57;
t5 = -t44 * t11 + t63 * t12;
t16 = -t30 * pkin(4) + qJD(5) - t19;
t47 = sin(qJ(6));
t43 = t49 ^ 2;
t42 = t48 ^ 2;
t40 = qJD(2) ^ 2 / 0.2e1;
t29 = qJD(6) + t32;
t25 = t63 * qJD(2) + t44 * t30;
t23 = t44 * qJD(2) - t63 * t30;
t18 = -qJD(2) * pkin(3) + t54;
t17 = t30 * pkin(3) + t52;
t15 = -t47 * t23 + t68 * t25;
t13 = t68 * t23 + t47 * t25;
t7 = t23 * pkin(5) + t16;
t4 = -t23 * pkin(8) + t6;
t3 = t32 * pkin(5) - t25 * pkin(8) + t5;
t2 = t47 * t3 + t68 * t4;
t1 = t68 * t3 - t47 * t4;
t8 = [0, 0, 0, 0, 0, t69, 0, 0, 0, 0, t42 * t69, t48 * t66, t56, t43 * t69, t55, t40, pkin(1) * t66 - pkin(7) * t56, -t50 * pkin(1) * t48 - pkin(7) * t55 (t42 + t43) * t50 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t43 / 0.2e1 + t42 / 0.2e1) * pkin(7) ^ 2) * t50, t28, -t67, t59, t58, -t60, t40, t20 * qJD(2) + t38 * t30, -t21 * qJD(2) + t38 * t32, -t20 * t32 - t21 * t30, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t40, -t59, t60, t28, -t67, t58, t18 * t32 + t19 * t30, t18 * qJD(2) - t17 * t30, -t19 * qJD(2) - t17 * t32, t17 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t32, t23 ^ 2 / 0.2e1, -t23 * t32, t28, t16 * t23 + t5 * t32, t16 * t25 - t6 * t32, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t29, t13 ^ 2 / 0.2e1, -t13 * t29, t29 ^ 2 / 0.2e1, t1 * t29 + t7 * t13, t7 * t15 - t2 * t29, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
