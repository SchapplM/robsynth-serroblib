% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:09
% EndTime: 2019-03-09 05:37:09
% DurationCPUTime: 0.14s
% Computational Cost: add. (395->54), mult. (758->115), div. (0->0), fcn. (415->6), ass. (0->46)
t43 = qJD(1) ^ 2;
t36 = t43 / 0.2e1;
t55 = -pkin(4) - pkin(5);
t54 = cos(qJ(4));
t53 = cos(qJ(6));
t40 = sin(qJ(4));
t42 = cos(qJ(3));
t49 = qJD(1) * t42;
t23 = -t54 * qJD(3) + t40 * t49;
t25 = t40 * qJD(3) + t54 * t49;
t52 = t23 * t25;
t41 = sin(qJ(3));
t32 = t41 * qJD(1) + qJD(4);
t51 = t23 * t32;
t19 = (pkin(3) * t41 - pkin(8) * t42 + qJ(2)) * qJD(1);
t31 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t20 = qJD(3) * pkin(8) + t41 * t31;
t10 = t40 * t19 + t54 * t20;
t50 = t43 * qJ(2);
t48 = qJD(3) * t31;
t47 = t23 ^ 2 / 0.2e1;
t46 = qJD(1) * qJD(3);
t7 = t32 * qJ(5) + t10;
t9 = t54 * t19 - t40 * t20;
t21 = -qJD(3) * pkin(3) - t42 * t31;
t45 = qJD(5) - t9;
t44 = t25 * qJ(5) - t21;
t39 = sin(qJ(6));
t38 = t42 ^ 2;
t37 = t41 ^ 2;
t34 = qJ(2) ^ 2 * t36;
t33 = -pkin(1) * qJD(1) + qJD(2);
t29 = -qJD(6) + t32;
t26 = t32 ^ 2 / 0.2e1;
t22 = t25 ^ 2 / 0.2e1;
t14 = t25 * t32;
t13 = t39 * t23 + t53 * t25;
t11 = -t53 * t23 + t39 * t25;
t8 = t23 * pkin(4) - t44;
t6 = -t32 * pkin(4) + t45;
t5 = t55 * t23 + t44;
t4 = t23 * pkin(9) + t7;
t3 = -t25 * pkin(9) + t55 * t32 + t45;
t2 = t39 * t3 + t53 * t4;
t1 = t53 * t3 - t39 * t4;
t12 = [0, 0, 0, 0, 0, t36, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, t33 * qJD(1), t50, t34 + t33 ^ 2 / 0.2e1, t38 * t36, -t42 * t43 * t41, t42 * t46, t37 * t36, -t41 * t46, qJD(3) ^ 2 / 0.2e1, t41 * t50 + t42 * t48, -t41 * t48 + t42 * t50 (-t37 - t38) * t31 * qJD(1), t34 + (t37 / 0.2e1 + t38 / 0.2e1) * t31 ^ 2, t22, -t52, t14, t47, -t51, t26, t21 * t23 + t9 * t32, -t10 * t32 + t21 * t25, -t10 * t23 - t9 * t25, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t22, t14, t52, t26, t51, t47, t8 * t23 - t6 * t32, -t7 * t23 + t6 * t25, -t8 * t25 + t7 * t32, t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, -t13 * t29, t11 ^ 2 / 0.2e1, t11 * t29, t29 ^ 2 / 0.2e1, -t1 * t29 + t5 * t11, t5 * t13 + t2 * t29, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
