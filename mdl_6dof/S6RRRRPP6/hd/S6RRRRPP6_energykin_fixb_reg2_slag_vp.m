% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPP6
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:15:03
% EndTime: 2019-03-09 21:15:03
% DurationCPUTime: 0.17s
% Computational Cost: add. (622->62), mult. (1356->127), div. (0->0), fcn. (898->6), ass. (0->48)
t49 = qJD(1) ^ 2;
t62 = t49 / 0.2e1;
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t28 = (-pkin(2) * t48 - pkin(8) * t46 - pkin(1)) * qJD(1);
t56 = t48 * qJD(1);
t35 = pkin(7) * t56 + qJD(2) * pkin(8);
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t22 = t47 * t28 - t45 * t35;
t57 = qJD(1) * t46;
t31 = t45 * qJD(2) + t47 * t57;
t38 = -qJD(3) + t56;
t10 = -t38 * pkin(3) - t31 * pkin(9) + t22;
t23 = t45 * t28 + t47 * t35;
t29 = -t47 * qJD(2) + t45 * t57;
t13 = -t29 * pkin(9) + t23;
t44 = sin(qJ(4));
t61 = cos(qJ(4));
t7 = t44 * t10 + t61 * t13;
t18 = t61 * t29 + t44 * t31;
t20 = -t44 * t29 + t61 * t31;
t60 = t18 * t20;
t36 = -qJD(4) + t38;
t14 = t18 * t36;
t15 = t20 * t36;
t59 = t48 * t49;
t58 = pkin(4) + qJ(6);
t16 = t18 ^ 2 / 0.2e1;
t17 = t20 ^ 2 / 0.2e1;
t55 = qJD(1) * qJD(2);
t54 = t48 * t55;
t53 = t46 * t55;
t34 = -qJD(2) * pkin(2) + pkin(7) * t57;
t5 = t36 * qJ(5) - t7;
t6 = t61 * t10 - t44 * t13;
t24 = t29 * pkin(3) + t34;
t52 = qJD(5) - t6;
t51 = -t20 * qJ(5) + t24;
t43 = t48 ^ 2;
t42 = t46 ^ 2;
t33 = t36 ^ 2 / 0.2e1;
t8 = t18 * pkin(4) + t51;
t4 = t36 * pkin(4) + t52;
t3 = t58 * t18 + t51;
t2 = -t18 * pkin(5) + qJD(6) - t5;
t1 = t20 * pkin(5) + t58 * t36 + t52;
t9 = [0, 0, 0, 0, 0, t62, 0, 0, 0, 0, t42 * t62, t46 * t59, t53, t43 * t62, t54, qJD(2) ^ 2 / 0.2e1, pkin(1) * t59 - pkin(7) * t53, -t49 * pkin(1) * t46 - pkin(7) * t54 (t42 + t43) * t49 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t43 / 0.2e1 + t42 / 0.2e1) * pkin(7) ^ 2) * t49, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t38, t29 ^ 2 / 0.2e1, t29 * t38, t38 ^ 2 / 0.2e1, -t22 * t38 + t34 * t29, t23 * t38 + t34 * t31, -t22 * t31 - t23 * t29, t23 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t17, -t60, -t15, t16, t14, t33, t24 * t18 - t6 * t36, t24 * t20 + t7 * t36, -t7 * t18 - t6 * t20, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t33, t15, -t14, t17, -t60, t16, t5 * t18 + t4 * t20, -t8 * t18 - t4 * t36, -t8 * t20 + t5 * t36, t8 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t33, -t14, -t15, t16, t60, t17, t1 * t20 - t2 * t18, -t2 * t36 - t3 * t20, t1 * t36 + t3 * t18, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg  = t9;
