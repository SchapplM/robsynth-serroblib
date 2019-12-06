% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:11
% EndTime: 2019-12-05 17:10:11
% DurationCPUTime: 0.14s
% Computational Cost: add. (207->31), mult. (350->91), div. (0->0), fcn. (209->8), ass. (0->36)
t25 = sin(qJ(4));
t22 = t25 ^ 2;
t41 = t22 / 0.2e1;
t28 = cos(qJ(4));
t23 = t28 ^ 2;
t40 = t23 / 0.2e1;
t39 = cos(qJ(5));
t21 = qJD(2) + qJD(3);
t38 = t21 * t25;
t37 = t21 * t28;
t30 = cos(qJ(2));
t16 = qJD(2) * pkin(2) + t30 * qJD(1);
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t27 = sin(qJ(2));
t36 = qJD(1) * t27;
t10 = t26 * t16 + t29 * t36;
t35 = qJD(4) * t25;
t34 = qJD(4) * t28;
t33 = qJD(1) * qJD(2);
t8 = t21 * pkin(7) + t10;
t32 = pkin(8) * t21 + t8;
t9 = t29 * t16 - t26 * t36;
t31 = qJD(1) ^ 2;
t24 = sin(qJ(5));
t20 = qJD(4) + qJD(5);
t19 = t21 ^ 2;
t13 = (t24 * t28 + t39 * t25) * t21;
t11 = t24 * t38 - t39 * t37;
t7 = -t21 * pkin(3) - t9;
t5 = (-pkin(4) * t28 - pkin(3)) * t21 - t9;
t4 = t32 * t28;
t3 = qJD(4) * pkin(4) - t32 * t25;
t2 = t24 * t3 + t39 * t4;
t1 = -t24 * t4 + t39 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t31 / 0.2e1, 0, 0, 0, 0, 0, qJD(2) ^ 2 / 0.2e1, t30 * t33, -t27 * t33, 0, (t27 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1) * t31, 0, 0, 0, 0, 0, t19 / 0.2e1, t9 * t21, -t10 * t21, 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t19 * t41, t25 * t19 * t28, t21 * t35, t19 * t40, t21 * t34, qJD(4) ^ 2 / 0.2e1, -t8 * t35 - t7 * t37, -t8 * t34 + t7 * t38, (t22 + t23) * t8 * t21, t7 ^ 2 / 0.2e1 + (t40 + t41) * t8 ^ 2, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t20, t11 ^ 2 / 0.2e1, -t11 * t20, t20 ^ 2 / 0.2e1, t1 * t20 + t5 * t11, t5 * t13 - t2 * t20, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
