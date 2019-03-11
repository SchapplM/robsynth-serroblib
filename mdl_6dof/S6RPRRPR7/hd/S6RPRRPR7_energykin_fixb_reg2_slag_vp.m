% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:20:59
% EndTime: 2019-03-09 05:21:00
% DurationCPUTime: 0.18s
% Computational Cost: add. (661->56), mult. (1342->131), div. (0->0), fcn. (870->8), ass. (0->44)
t50 = qJD(1) ^ 2;
t40 = t50 / 0.2e1;
t49 = cos(qJ(3));
t33 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t51 = -pkin(8) * qJD(1) + t33;
t25 = qJD(3) * pkin(3) + t51 * t49;
t47 = sin(qJ(3));
t27 = t51 * t47;
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t16 = t46 * t25 + t48 * t27;
t28 = (-t46 * t49 - t47 * t48) * qJD(1);
t11 = t28 * qJ(5) + t16;
t44 = sin(pkin(10));
t54 = cos(pkin(10));
t15 = t48 * t25 - t46 * t27;
t30 = (-t46 * t47 + t48 * t49) * qJD(1);
t39 = qJD(3) + qJD(4);
t9 = t39 * pkin(4) - t30 * qJ(5) + t15;
t6 = t54 * t11 + t44 * t9;
t56 = cos(qJ(6));
t31 = (pkin(3) * t47 + qJ(2)) * qJD(1);
t55 = t50 * qJ(2);
t53 = qJD(3) * t33;
t52 = qJD(1) * qJD(3);
t18 = -t54 * t28 + t44 * t30;
t21 = -t28 * pkin(4) + qJD(5) + t31;
t5 = -t44 * t11 + t54 * t9;
t45 = sin(qJ(6));
t43 = t49 ^ 2;
t42 = t47 ^ 2;
t38 = qJ(2) ^ 2 * t40;
t37 = -pkin(1) * qJD(1) + qJD(2);
t35 = t39 ^ 2 / 0.2e1;
t20 = t44 * t28 + t54 * t30;
t17 = qJD(6) + t18;
t14 = t56 * t20 + t45 * t39;
t12 = t45 * t20 - t56 * t39;
t7 = t18 * pkin(5) - t20 * pkin(9) + t21;
t4 = t39 * pkin(9) + t6;
t3 = -t39 * pkin(5) - t5;
t2 = t56 * t4 + t45 * t7;
t1 = -t45 * t4 + t56 * t7;
t8 = [0, 0, 0, 0, 0, t40, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, t37 * qJD(1), t55, t38 + t37 ^ 2 / 0.2e1, t43 * t40, -t49 * t50 * t47, t49 * t52, t42 * t40, -t47 * t52, qJD(3) ^ 2 / 0.2e1, t47 * t55 + t49 * t53, -t47 * t53 + t49 * t55 (-t42 - t43) * t33 * qJD(1), t38 + (t42 / 0.2e1 + t43 / 0.2e1) * t33 ^ 2, t30 ^ 2 / 0.2e1, t30 * t28, t30 * t39, t28 ^ 2 / 0.2e1, t28 * t39, t35, t15 * t39 - t31 * t28, -t16 * t39 + t31 * t30, -t15 * t30 + t16 * t28, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t39, t18 ^ 2 / 0.2e1, -t18 * t39, t35, t21 * t18 + t5 * t39, t21 * t20 - t6 * t39, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t17, t12 ^ 2 / 0.2e1, -t12 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t12, t3 * t14 - t2 * t17, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
