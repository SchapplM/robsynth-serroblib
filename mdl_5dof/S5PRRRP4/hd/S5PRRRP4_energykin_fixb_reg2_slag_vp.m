% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRP4
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:19
% EndTime: 2019-12-05 16:46:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (138->28), mult. (250->73), div. (0->0), fcn. (128->6), ass. (0->34)
t21 = sin(qJ(4));
t19 = t21 ^ 2;
t37 = t19 / 0.2e1;
t24 = cos(qJ(4));
t20 = t24 ^ 2;
t36 = t20 / 0.2e1;
t26 = cos(qJ(2));
t10 = qJD(2) * pkin(2) + t26 * qJD(1);
t22 = sin(qJ(3));
t25 = cos(qJ(3));
t23 = sin(qJ(2));
t33 = qJD(1) * t23;
t8 = t22 * t10 + t25 * t33;
t17 = qJD(2) + qJD(3);
t35 = t17 * t21;
t34 = t17 * t24;
t32 = qJD(4) * t21;
t31 = qJD(4) * t24;
t30 = qJD(1) * qJD(2);
t16 = t17 ^ 2;
t29 = t21 * t16 * t24;
t28 = t17 * t31;
t7 = t25 * t10 - t22 * t33;
t27 = qJD(1) ^ 2;
t18 = qJD(4) ^ 2 / 0.2e1;
t13 = t17 * t32;
t12 = t16 * t36;
t11 = t16 * t37;
t6 = t17 * pkin(7) + t8;
t5 = -t17 * pkin(3) - t7;
t3 = qJD(4) * qJ(5) + t24 * t6;
t2 = -qJD(4) * pkin(4) + t21 * t6 + qJD(5);
t1 = (-pkin(4) * t24 - qJ(5) * t21 - pkin(3)) * t17 - t7;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t27 / 0.2e1, 0, 0, 0, 0, 0, qJD(2) ^ 2 / 0.2e1, t26 * t30, -t23 * t30, 0, (t23 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1) * t27, 0, 0, 0, 0, 0, t16 / 0.2e1, t7 * t17, -t8 * t17, 0, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t11, t29, t13, t12, t28, t18, -t6 * t32 - t5 * t34, -t6 * t31 + t5 * t35, (t19 + t20) * t6 * t17, t5 ^ 2 / 0.2e1 + (t36 + t37) * t6 ^ 2, t11, t13, -t29, t18, -t28, t12, -t2 * qJD(4) - t1 * t34, (t2 * t21 + t24 * t3) * t17, t3 * qJD(4) - t1 * t35, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t4;
