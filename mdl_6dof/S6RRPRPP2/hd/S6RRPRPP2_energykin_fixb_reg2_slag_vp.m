% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:27
% EndTime: 2019-03-09 09:52:28
% DurationCPUTime: 0.18s
% Computational Cost: add. (533->59), mult. (1310->123), div. (0->0), fcn. (879->6), ass. (0->49)
t47 = qJD(1) ^ 2;
t63 = t47 / 0.2e1;
t62 = -pkin(4) - pkin(5);
t61 = cos(qJ(4));
t43 = sin(pkin(9));
t45 = sin(qJ(2));
t46 = cos(qJ(2));
t57 = cos(pkin(9));
t32 = (t43 * t46 + t57 * t45) * qJD(1);
t44 = sin(qJ(4));
t21 = -t61 * qJD(2) + t44 * t32;
t23 = t44 * qJD(2) + t61 * t32;
t9 = t23 * t21;
t55 = qJD(1) * t46;
t56 = qJD(1) * t45;
t30 = t43 * t56 - t57 * t55;
t28 = qJD(4) + t30;
t12 = t23 * t28;
t60 = t28 * t21;
t59 = t46 * t47;
t58 = pkin(7) + qJ(3);
t37 = qJD(3) + (-pkin(2) * t46 - pkin(1)) * qJD(1);
t11 = t30 * pkin(3) - t32 * pkin(8) + t37;
t35 = qJD(2) * pkin(2) - t58 * t56;
t36 = t58 * t55;
t18 = t43 * t35 + t57 * t36;
t16 = qJD(2) * pkin(8) + t18;
t7 = t44 * t11 + t61 * t16;
t17 = t57 * t35 - t43 * t36;
t19 = t21 ^ 2 / 0.2e1;
t24 = t28 ^ 2 / 0.2e1;
t54 = qJD(1) * qJD(2);
t5 = t28 * qJ(5) + t7;
t53 = t45 * t54;
t52 = t46 * t54;
t51 = qJD(2) * pkin(3) + t17;
t6 = t61 * t11 - t44 * t16;
t50 = qJD(5) - t6;
t49 = t23 * qJ(5) + t51;
t42 = t46 ^ 2;
t41 = t45 ^ 2;
t40 = qJD(2) ^ 2 / 0.2e1;
t20 = t23 ^ 2 / 0.2e1;
t8 = t21 * pkin(4) - t49;
t4 = -t28 * pkin(4) + t50;
t3 = t62 * t21 + qJD(6) + t49;
t2 = t21 * qJ(6) + t5;
t1 = -t23 * qJ(6) + t62 * t28 + t50;
t10 = [0, 0, 0, 0, 0, t63, 0, 0, 0, 0, t41 * t63, t45 * t59, t53, t42 * t63, t52, t40, pkin(1) * t59 - pkin(7) * t53, -t47 * pkin(1) * t45 - pkin(7) * t52 (t41 + t42) * t47 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * pkin(7) ^ 2) * t47, t32 ^ 2 / 0.2e1, -t32 * t30, t32 * qJD(2), t30 ^ 2 / 0.2e1, -t30 * qJD(2), t40, t17 * qJD(2) + t37 * t30, -t18 * qJD(2) + t37 * t32, -t17 * t32 - t18 * t30, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t20, -t9, t12, t19, -t60, t24, -t21 * t51 + t6 * t28, -t23 * t51 - t7 * t28, -t7 * t21 - t6 * t23, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1, t20, t12, t9, t24, t60, t19, t8 * t21 - t4 * t28, -t5 * t21 + t4 * t23, -t8 * t23 + t5 * t28, t5 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t20, t9, -t12, t19, -t60, t24, -t1 * t28 - t3 * t21, t2 * t28 + t3 * t23, -t1 * t23 + t2 * t21, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t10;
