% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:09:06
% EndTime: 2019-03-09 21:09:06
% DurationCPUTime: 0.19s
% Computational Cost: add. (622->60), mult. (1356->127), div. (0->0), fcn. (898->6), ass. (0->48)
t49 = qJD(1) ^ 2;
t63 = t49 / 0.2e1;
t62 = pkin(4) + pkin(5);
t61 = cos(qJ(4));
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t46 = sin(qJ(2));
t58 = qJD(1) * t46;
t28 = -t47 * qJD(2) + t45 * t58;
t30 = t45 * qJD(2) + t47 * t58;
t44 = sin(qJ(4));
t18 = t61 * t28 + t44 * t30;
t48 = cos(qJ(2));
t57 = t48 * qJD(1);
t38 = -qJD(3) + t57;
t36 = -qJD(4) + t38;
t60 = t18 * t36;
t20 = -t44 * t28 + t61 * t30;
t9 = t20 * t18;
t15 = t20 * t36;
t59 = t48 * t49;
t27 = (-pkin(2) * t48 - pkin(8) * t46 - pkin(1)) * qJD(1);
t35 = pkin(7) * t57 + qJD(2) * pkin(8);
t21 = t47 * t27 - t45 * t35;
t11 = -t38 * pkin(3) - t30 * pkin(9) + t21;
t22 = t45 * t27 + t47 * t35;
t14 = -t28 * pkin(9) + t22;
t7 = t44 * t11 + t61 * t14;
t16 = t18 ^ 2 / 0.2e1;
t32 = t36 ^ 2 / 0.2e1;
t56 = qJD(1) * qJD(2);
t5 = -t36 * qJ(5) + t7;
t55 = t46 * t56;
t54 = t48 * t56;
t34 = -qJD(2) * pkin(2) + pkin(7) * t58;
t6 = t61 * t11 - t44 * t14;
t53 = -t28 * pkin(3) - t34;
t52 = qJD(5) - t6;
t51 = t20 * qJ(5) + t53;
t43 = t48 ^ 2;
t42 = t46 ^ 2;
t17 = t20 ^ 2 / 0.2e1;
t8 = t18 * pkin(4) - t51;
t4 = t36 * pkin(4) + t52;
t3 = -t62 * t18 + qJD(6) + t51;
t2 = t18 * qJ(6) + t5;
t1 = -t20 * qJ(6) + t62 * t36 + t52;
t10 = [0, 0, 0, 0, 0, t63, 0, 0, 0, 0, t42 * t63, t46 * t59, t55, t43 * t63, t54, qJD(2) ^ 2 / 0.2e1, pkin(1) * t59 - pkin(7) * t55, -t49 * pkin(1) * t46 - pkin(7) * t54 (t42 + t43) * t49 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t43 / 0.2e1 + t42 / 0.2e1) * pkin(7) ^ 2) * t49, t30 ^ 2 / 0.2e1, -t30 * t28, -t30 * t38, t28 ^ 2 / 0.2e1, t28 * t38, t38 ^ 2 / 0.2e1, -t21 * t38 + t34 * t28, t22 * t38 + t34 * t30, -t21 * t30 - t22 * t28, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t17, -t9, -t15, t16, t60, t32, -t18 * t53 - t6 * t36, -t20 * t53 + t7 * t36, -t7 * t18 - t6 * t20, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t53 ^ 2 / 0.2e1, t17, -t15, t9, t32, -t60, t16, t8 * t18 + t4 * t36, -t5 * t18 + t4 * t20, -t8 * t20 - t5 * t36, t5 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t17, t9, t15, t16, t60, t32, t1 * t36 - t3 * t18, -t2 * t36 + t3 * t20, -t1 * t20 + t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t10;
