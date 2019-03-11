% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:09
% EndTime: 2019-03-09 01:56:09
% DurationCPUTime: 0.16s
% Computational Cost: add. (356->49), mult. (766->109), div. (0->0), fcn. (452->6), ass. (0->44)
t41 = qJD(1) ^ 2;
t35 = t41 / 0.2e1;
t54 = pkin(4) + pkin(8);
t53 = cos(qJ(6));
t52 = sin(qJ(4));
t37 = sin(pkin(9));
t38 = cos(pkin(9));
t40 = cos(qJ(4));
t21 = (-t37 * t40 - t52 * t38) * qJD(1);
t49 = qJD(1) * t38;
t50 = qJD(1) * t37;
t23 = t40 * t49 - t52 * t50;
t51 = t23 * t21;
t27 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t44 = -pkin(7) * qJD(1) + t27;
t18 = t44 * t37;
t19 = t44 * t38;
t10 = t40 * t18 + t52 * t19;
t48 = qJD(4) * t21;
t47 = t23 * qJD(4);
t46 = t21 ^ 2 / 0.2e1;
t45 = t23 ^ 2 / 0.2e1;
t30 = qJD(1) * qJ(2) + qJD(3);
t25 = pkin(3) * t50 + t30;
t9 = -t52 * t18 + t40 * t19;
t43 = qJD(5) - t9;
t7 = -qJD(4) * qJ(5) - t10;
t42 = -t23 * qJ(5) + t25;
t39 = sin(qJ(6));
t34 = qJD(4) ^ 2 / 0.2e1;
t33 = t38 ^ 2;
t32 = t37 ^ 2;
t31 = -pkin(1) * qJD(1) + qJD(2);
t20 = qJD(6) + t23;
t13 = t53 * qJD(4) - t39 * t21;
t11 = t39 * qJD(4) + t53 * t21;
t8 = -t21 * pkin(4) + t42;
t6 = -qJD(4) * pkin(4) + t43;
t5 = -t54 * t21 + t42;
t4 = t21 * pkin(5) - t7;
t3 = t23 * pkin(5) - t54 * qJD(4) + t43;
t2 = t39 * t3 + t53 * t5;
t1 = t53 * t3 - t39 * t5;
t12 = [0, 0, 0, 0, 0, t35, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, t31 * qJD(1), t41 * qJ(2), qJ(2) ^ 2 * t35 + t31 ^ 2 / 0.2e1, t33 * t35, -t38 * t41 * t37, 0, t32 * t35, 0, 0, t30 * t50, t30 * t49 (-t32 - t33) * t27 * qJD(1), t30 ^ 2 / 0.2e1 + (t32 / 0.2e1 + t33 / 0.2e1) * t27 ^ 2, t45, t51, t47, t46, t48, t34, t9 * qJD(4) - t25 * t21, -t10 * qJD(4) + t25 * t23, t10 * t21 - t9 * t23, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t34, -t47, -t48, t45, t51, t46, -t7 * t21 + t6 * t23, t6 * qJD(4) + t8 * t21, -t7 * qJD(4) - t8 * t23, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t20, t11 ^ 2 / 0.2e1, -t11 * t20, t20 ^ 2 / 0.2e1, t1 * t20 + t4 * t11, t4 * t13 - t2 * t20, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t12;
