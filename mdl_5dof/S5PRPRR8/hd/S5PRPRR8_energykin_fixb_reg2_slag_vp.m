% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:25
% EndTime: 2019-12-05 16:04:25
% DurationCPUTime: 0.15s
% Computational Cost: add. (135->36), mult. (316->89), div. (0->0), fcn. (180->8), ass. (0->37)
t29 = qJD(2) ^ 2;
t21 = t29 / 0.2e1;
t30 = qJD(1) ^ 2;
t42 = t30 / 0.2e1;
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t23 = cos(pkin(5));
t39 = qJD(1) * t23;
t28 = cos(qJ(2));
t22 = sin(pkin(5));
t40 = qJD(1) * t22;
t31 = -t28 * t40 + qJD(3);
t9 = (-pkin(2) - pkin(7)) * qJD(2) + t31;
t6 = t25 * t9 + t27 * t39;
t41 = cos(qJ(5));
t38 = qJD(2) * t27;
t26 = sin(qJ(2));
t33 = t26 * t40;
t14 = qJD(2) * qJ(3) + t33;
t37 = t14 * qJD(2);
t36 = t25 * qJD(2);
t35 = t14 ^ 2 / 0.2e1;
t34 = qJD(2) * qJD(4);
t32 = qJD(2) * t40;
t5 = -t25 * t39 + t27 * t9;
t24 = sin(qJ(5));
t18 = t23 ^ 2 * t42;
t17 = qJD(5) + t36;
t13 = t24 * qJD(4) + t41 * t38;
t11 = -t41 * qJD(4) + t24 * t38;
t10 = -qJD(2) * pkin(2) + t31;
t7 = t33 + (pkin(4) * t25 - pkin(8) * t27 + qJ(3)) * qJD(2);
t4 = qJD(4) * pkin(8) + t6;
t3 = -qJD(4) * pkin(4) - t5;
t2 = t24 * t7 + t41 * t4;
t1 = -t24 * t4 + t41 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, t21, t28 * t32, -t26 * t32, 0, t18 + (t26 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1) * t30 * t22 ^ 2, t21, 0, 0, 0, 0, 0, 0, t10 * qJD(2), t37, t18 + t35 + t10 ^ 2 / 0.2e1, t27 ^ 2 * t21, -t27 * t29 * t25, t27 * t34, t25 ^ 2 * t21, -t25 * t34, qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) + t14 * t36, -t6 * qJD(4) + t27 * t37, (-t25 * t6 - t27 * t5) * qJD(2), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t35, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t17, t11 ^ 2 / 0.2e1, -t11 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t11, t3 * t13 - t2 * t17, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
