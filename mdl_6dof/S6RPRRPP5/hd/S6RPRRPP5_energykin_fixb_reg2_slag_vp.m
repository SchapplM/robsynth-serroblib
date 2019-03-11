% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:40
% EndTime: 2019-03-09 04:44:40
% DurationCPUTime: 0.18s
% Computational Cost: add. (502->58), mult. (1264->118), div. (0->0), fcn. (875->6), ass. (0->45)
t47 = qJD(1) ^ 2;
t59 = t47 / 0.2e1;
t58 = -pkin(4) - pkin(5);
t57 = cos(qJ(3));
t56 = cos(qJ(4));
t43 = sin(pkin(9));
t44 = cos(pkin(9));
t46 = sin(qJ(3));
t34 = (t57 * t43 + t44 * t46) * qJD(1);
t45 = sin(qJ(4));
t21 = -t56 * qJD(3) + t45 * t34;
t23 = t45 * qJD(3) + t56 * t34;
t9 = t23 * t21;
t52 = qJD(1) * t44;
t53 = qJD(1) * t43;
t32 = t46 * t53 - t57 * t52;
t26 = qJD(4) + t32;
t12 = t23 * t26;
t55 = t26 * t21;
t54 = pkin(7) + qJ(2);
t37 = qJD(2) + (-pkin(2) * t44 - pkin(1)) * qJD(1);
t11 = t32 * pkin(3) - t34 * pkin(8) + t37;
t35 = t54 * t53;
t36 = t54 * t52;
t18 = -t46 * t35 + t57 * t36;
t16 = qJD(3) * pkin(8) + t18;
t7 = t45 * t11 + t56 * t16;
t17 = -t57 * t35 - t46 * t36;
t19 = t21 ^ 2 / 0.2e1;
t24 = t26 ^ 2 / 0.2e1;
t5 = t26 * qJ(5) + t7;
t51 = qJD(3) * pkin(3) + t17;
t6 = t56 * t11 - t45 * t16;
t50 = qJD(5) - t6;
t49 = t23 * qJ(5) + t51;
t42 = t44 ^ 2;
t41 = t43 ^ 2;
t39 = -qJD(1) * pkin(1) + qJD(2);
t20 = t23 ^ 2 / 0.2e1;
t8 = t21 * pkin(4) - t49;
t4 = -t26 * pkin(4) + t50;
t3 = t58 * t21 + qJD(6) + t49;
t2 = t21 * qJ(6) + t5;
t1 = -t23 * qJ(6) + t58 * t26 + t50;
t10 = [0, 0, 0, 0, 0, t59, 0, 0, 0, 0, t41 * t59, t43 * t47 * t44, 0, t42 * t59, 0, 0, -t39 * t52, t39 * t53 (t41 + t42) * t47 * qJ(2), t39 ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * qJ(2) ^ 2 * t47, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * qJD(3), t32 ^ 2 / 0.2e1, -t32 * qJD(3), qJD(3) ^ 2 / 0.2e1, t17 * qJD(3) + t37 * t32, -t18 * qJD(3) + t37 * t34, -t17 * t34 - t18 * t32, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t20, -t9, t12, t19, -t55, t24, -t21 * t51 + t6 * t26, -t23 * t51 - t7 * t26, -t7 * t21 - t6 * t23, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1, t20, t12, t9, t24, t55, t19, t8 * t21 - t4 * t26, -t5 * t21 + t4 * t23, -t8 * t23 + t5 * t26, t5 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t20, t9, -t12, t19, -t55, t24, -t1 * t26 - t3 * t21, t2 * t26 + t3 * t23, -t1 * t23 + t2 * t21, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t10;
