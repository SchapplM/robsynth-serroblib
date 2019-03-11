% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:02
% EndTime: 2019-03-09 03:56:02
% DurationCPUTime: 0.18s
% Computational Cost: add. (644->56), mult. (1342->131), div. (0->0), fcn. (870->8), ass. (0->44)
t50 = qJD(1) ^ 2;
t40 = t50 / 0.2e1;
t49 = cos(qJ(3));
t33 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t51 = -qJ(4) * qJD(1) + t33;
t25 = qJD(3) * pkin(3) + t51 * t49;
t48 = sin(qJ(3));
t27 = t51 * t48;
t44 = sin(pkin(10));
t45 = cos(pkin(10));
t15 = t45 * t25 - t44 * t27;
t30 = (-t44 * t48 + t45 * t49) * qJD(1);
t10 = qJD(3) * pkin(4) - t30 * pkin(8) + t15;
t16 = t44 * t25 + t45 * t27;
t28 = (-t44 * t49 - t45 * t48) * qJD(1);
t11 = t28 * pkin(8) + t16;
t47 = sin(qJ(5));
t56 = cos(qJ(5));
t6 = t47 * t10 + t56 * t11;
t55 = cos(qJ(6));
t54 = t50 * qJ(2);
t53 = qJD(3) * t33;
t52 = qJD(1) * qJD(3);
t31 = qJD(4) + (pkin(3) * t48 + qJ(2)) * qJD(1);
t18 = -t56 * t28 + t47 * t30;
t21 = -t28 * pkin(4) + t31;
t5 = t56 * t10 - t47 * t11;
t46 = sin(qJ(6));
t43 = t49 ^ 2;
t42 = t48 ^ 2;
t39 = qJD(3) ^ 2 / 0.2e1;
t38 = qJD(3) + qJD(5);
t37 = qJ(2) ^ 2 * t40;
t36 = -pkin(1) * qJD(1) + qJD(2);
t20 = t47 * t28 + t56 * t30;
t17 = qJD(6) + t18;
t14 = t55 * t20 + t46 * t38;
t12 = t46 * t20 - t55 * t38;
t7 = t18 * pkin(5) - t20 * pkin(9) + t21;
t4 = t38 * pkin(9) + t6;
t3 = -t38 * pkin(5) - t5;
t2 = t55 * t4 + t46 * t7;
t1 = -t46 * t4 + t55 * t7;
t8 = [0, 0, 0, 0, 0, t40, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, t36 * qJD(1), t54, t37 + t36 ^ 2 / 0.2e1, t43 * t40, -t49 * t50 * t48, t49 * t52, t42 * t40, -t48 * t52, t39, t48 * t54 + t49 * t53, -t48 * t53 + t49 * t54 (-t42 - t43) * t33 * qJD(1), t37 + (t42 / 0.2e1 + t43 / 0.2e1) * t33 ^ 2, t30 ^ 2 / 0.2e1, t30 * t28, t30 * qJD(3), t28 ^ 2 / 0.2e1, t28 * qJD(3), t39, t15 * qJD(3) - t31 * t28, -t16 * qJD(3) + t31 * t30, -t15 * t30 + t16 * t28, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t38, t18 ^ 2 / 0.2e1, -t18 * t38, t38 ^ 2 / 0.2e1, t21 * t18 + t5 * t38, t21 * t20 - t6 * t38, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t17, t12 ^ 2 / 0.2e1, -t12 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t12, t3 * t14 - t2 * t17, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
