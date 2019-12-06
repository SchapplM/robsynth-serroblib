% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRP3
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:12
% EndTime: 2019-12-05 16:44:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (156->36), mult. (405->82), div. (0->0), fcn. (240->4), ass. (0->32)
t32 = qJD(2) ^ 2;
t38 = t32 / 0.2e1;
t31 = cos(qJ(3));
t26 = t31 * qJD(1);
t30 = sin(qJ(3));
t35 = qJD(2) * t30;
t10 = qJD(3) * pkin(3) + t26 + (-pkin(7) - pkin(6)) * t35;
t34 = qJD(2) * t31;
t20 = pkin(6) * t34 + t30 * qJD(1);
t15 = pkin(7) * t34 + t20;
t29 = sin(qJ(4));
t37 = cos(qJ(4));
t4 = t29 * t10 + t37 * t15;
t36 = t31 * t32;
t33 = qJD(2) * qJD(3);
t3 = t37 * t10 - t29 * t15;
t21 = (-pkin(3) * t31 - pkin(2)) * qJD(2);
t28 = qJD(1) ^ 2 / 0.2e1;
t27 = qJD(3) + qJD(4);
t23 = t27 ^ 2 / 0.2e1;
t19 = -pkin(6) * t35 + t26;
t18 = (t29 * t31 + t37 * t30) * qJD(2);
t16 = t29 * t35 - t37 * t34;
t14 = t18 ^ 2 / 0.2e1;
t13 = t16 ^ 2 / 0.2e1;
t12 = t18 * t27;
t11 = t16 * t27;
t6 = t16 * pkin(4) + qJD(5) + t21;
t5 = t18 * t16;
t2 = -t16 * qJ(5) + t4;
t1 = t27 * pkin(4) - t18 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, t38, 0, 0, 0, t28, t30 ^ 2 * t38, t30 * t36, t30 * t33, t31 ^ 2 * t38, t31 * t33, qJD(3) ^ 2 / 0.2e1, pkin(2) * t36 + t19 * qJD(3), -t32 * pkin(2) * t30 - t20 * qJD(3), (-t19 * t30 + t20 * t31) * qJD(2), t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + pkin(2) ^ 2 * t38, t14, -t5, t12, t13, -t11, t23, t21 * t16 + t3 * t27, t21 * t18 - t4 * t27, -t4 * t16 - t3 * t18, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14, -t5, t12, t13, -t11, t23, t1 * t27 + t6 * t16, t6 * t18 - t2 * t27, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg = t7;
