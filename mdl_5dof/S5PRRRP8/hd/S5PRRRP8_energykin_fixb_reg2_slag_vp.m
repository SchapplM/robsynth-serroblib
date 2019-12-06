% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:34
% EndTime: 2019-12-05 17:00:34
% DurationCPUTime: 0.17s
% Computational Cost: add. (192->41), mult. (486->99), div. (0->0), fcn. (311->8), ass. (0->40)
t35 = qJD(2) ^ 2;
t48 = t35 / 0.2e1;
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t34 = cos(qJ(2));
t28 = sin(pkin(5));
t44 = qJD(1) * t28;
t38 = t34 * t44;
t12 = -t38 + (-pkin(3) * t33 - pkin(8) * t31 - pkin(2)) * qJD(2);
t30 = sin(qJ(4));
t47 = cos(qJ(4));
t32 = sin(qJ(2));
t20 = qJD(2) * pkin(7) + t32 * t44;
t29 = cos(pkin(5));
t43 = qJD(1) * t29;
t11 = t33 * t20 + t31 * t43;
t8 = qJD(3) * pkin(8) + t11;
t4 = t30 * t12 + t47 * t8;
t42 = qJD(2) * t31;
t17 = -qJD(3) * t47 + t30 * t42;
t19 = t30 * qJD(3) + t42 * t47;
t46 = t19 * t17;
t41 = t33 * qJD(2);
t24 = -qJD(4) + t41;
t45 = t24 * t17;
t40 = t17 ^ 2 / 0.2e1;
t39 = qJD(2) * qJD(3);
t37 = qJD(2) * t44;
t10 = -t31 * t20 + t33 * t43;
t3 = t12 * t47 - t30 * t8;
t7 = -qJD(3) * pkin(3) - t10;
t36 = qJD(1) ^ 2;
t22 = t24 ^ 2 / 0.2e1;
t21 = -qJD(2) * pkin(2) - t38;
t14 = t19 ^ 2 / 0.2e1;
t13 = t19 * t24;
t5 = pkin(4) * t17 - qJ(5) * t19 + t7;
t2 = -qJ(5) * t24 + t4;
t1 = t24 * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t36 / 0.2e1, 0, 0, 0, 0, 0, t48, t34 * t37, -t32 * t37, 0, (t29 ^ 2 / 0.2e1 + (t32 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1) * t28 ^ 2) * t36, t31 ^ 2 * t48, t31 * t35 * t33, t31 * t39, t33 ^ 2 * t48, t33 * t39, qJD(3) ^ 2 / 0.2e1, qJD(3) * t10 - t21 * t41, -qJD(3) * t11 + t21 * t42, (-t10 * t31 + t11 * t33) * qJD(2), t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14, -t46, -t13, t40, t45, t22, t17 * t7 - t24 * t3, t19 * t7 + t24 * t4, -t17 * t4 - t19 * t3, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t14, -t13, t46, t22, -t45, t40, t1 * t24 + t17 * t5, t1 * t19 - t17 * t2, -t19 * t5 - t2 * t24, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
