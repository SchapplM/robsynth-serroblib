% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRP7
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:24
% EndTime: 2019-12-05 16:56:25
% DurationCPUTime: 0.16s
% Computational Cost: add. (192->41), mult. (492->99), div. (0->0), fcn. (317->8), ass. (0->40)
t39 = qJD(2) ^ 2;
t49 = t39 / 0.2e1;
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t38 = cos(qJ(2));
t32 = sin(pkin(5));
t47 = qJD(1) * t32;
t42 = t38 * t47;
t14 = -t42 + (-pkin(3) * t37 - pkin(8) * t35 - pkin(2)) * qJD(2);
t34 = sin(qJ(4));
t48 = cos(qJ(4));
t36 = sin(qJ(2));
t24 = qJD(2) * pkin(7) + t36 * t47;
t33 = cos(pkin(5));
t46 = qJD(1) * t33;
t13 = t37 * t24 + t35 * t46;
t9 = qJD(3) * pkin(8) + t13;
t4 = t34 * t14 + t48 * t9;
t45 = qJD(2) * t35;
t44 = t37 * qJD(2);
t43 = qJD(2) * qJD(3);
t3 = t48 * t14 - t34 * t9;
t41 = qJD(2) * t47;
t12 = -t35 * t24 + t37 * t46;
t8 = -qJD(3) * pkin(3) - t12;
t40 = qJD(1) ^ 2;
t28 = -qJD(4) + t44;
t26 = t28 ^ 2 / 0.2e1;
t25 = -qJD(2) * pkin(2) - t42;
t23 = t34 * qJD(3) + t48 * t45;
t21 = -t48 * qJD(3) + t34 * t45;
t18 = t23 ^ 2 / 0.2e1;
t17 = t21 ^ 2 / 0.2e1;
t16 = t23 * t28;
t15 = t21 * t28;
t7 = t23 * t21;
t5 = t21 * pkin(4) + qJD(5) + t8;
t2 = -t21 * qJ(5) + t4;
t1 = -t28 * pkin(4) - t23 * qJ(5) + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t40 / 0.2e1, 0, 0, 0, 0, 0, t49, t38 * t41, -t36 * t41, 0, (t33 ^ 2 / 0.2e1 + (t36 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1) * t32 ^ 2) * t40, t35 ^ 2 * t49, t35 * t39 * t37, t35 * t43, t37 ^ 2 * t49, t37 * t43, qJD(3) ^ 2 / 0.2e1, t12 * qJD(3) - t25 * t44, -t13 * qJD(3) + t25 * t45, (-t12 * t35 + t13 * t37) * qJD(2), t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t18, -t7, -t16, t17, t15, t26, t8 * t21 - t3 * t28, t8 * t23 + t4 * t28, -t4 * t21 - t3 * t23, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t18, -t7, -t16, t17, t15, t26, -t1 * t28 + t5 * t21, t2 * t28 + t5 * t23, -t1 * t23 - t2 * t21, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
