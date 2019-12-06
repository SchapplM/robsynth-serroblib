% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:16
% EndTime: 2019-12-05 16:52:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (167->35), mult. (413->91), div. (0->0), fcn. (239->6), ass. (0->38)
t29 = qJD(2) ^ 2;
t42 = t29 / 0.2e1;
t27 = cos(qJ(3));
t26 = sin(qJ(2));
t17 = qJD(2) * pkin(6) + t26 * qJD(1);
t31 = pkin(7) * qJD(2) + t17;
t10 = t31 * t27;
t24 = sin(qJ(4));
t41 = cos(qJ(4));
t25 = sin(qJ(3));
t8 = qJD(3) * pkin(3) - t31 * t25;
t5 = t41 * t10 + t24 * t8;
t37 = qJD(2) * t27;
t38 = qJD(2) * t25;
t12 = t24 * t38 - t41 * t37;
t14 = (t24 * t27 + t41 * t25) * qJD(2);
t40 = t14 * t12;
t21 = qJD(3) + qJD(4);
t39 = t21 * t12;
t36 = qJD(3) * t17;
t28 = cos(qJ(2));
t35 = t28 * qJD(1);
t34 = t12 ^ 2 / 0.2e1;
t33 = qJD(1) * qJD(2);
t32 = qJD(2) * qJD(3);
t4 = -t24 * t10 + t41 * t8;
t15 = -t35 + (-pkin(3) * t27 - pkin(2)) * qJD(2);
t30 = qJD(1) ^ 2;
t23 = t27 ^ 2;
t22 = t25 ^ 2;
t20 = t21 ^ 2 / 0.2e1;
t18 = -qJD(2) * pkin(2) - t35;
t11 = t14 ^ 2 / 0.2e1;
t9 = t14 * t21;
t3 = t12 * pkin(4) - t14 * qJ(5) + t15;
t2 = t21 * qJ(5) + t5;
t1 = -t21 * pkin(4) + qJD(5) - t4;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t30 / 0.2e1, 0, 0, 0, 0, 0, t42, t28 * t33, -t26 * t33, 0, (t26 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1) * t30, t22 * t42, t25 * t29 * t27, t25 * t32, t23 * t42, t27 * t32, qJD(3) ^ 2 / 0.2e1, -t18 * t37 - t25 * t36, t18 * t38 - t27 * t36, (t22 + t23) * t17 * qJD(2), t18 ^ 2 / 0.2e1 + (t23 / 0.2e1 + t22 / 0.2e1) * t17 ^ 2, t11, -t40, t9, t34, -t39, t20, t15 * t12 + t4 * t21, t15 * t14 - t5 * t21, -t5 * t12 - t4 * t14, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t11, t9, t40, t20, t39, t34, -t1 * t21 + t3 * t12, t1 * t14 - t2 * t12, -t3 * t14 + t2 * t21, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
