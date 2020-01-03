% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:55
% EndTime: 2019-12-31 19:52:55
% DurationCPUTime: 0.10s
% Computational Cost: add. (143->27), mult. (206->66), div. (0->0), fcn. (62->4), ass. (0->34)
t18 = sin(qJ(4));
t16 = t18 ^ 2;
t36 = t16 / 0.2e1;
t20 = cos(qJ(4));
t17 = t20 ^ 2;
t35 = t17 / 0.2e1;
t14 = qJD(1) + qJD(2);
t34 = t14 * t18;
t33 = t14 * t20;
t32 = pkin(1) * qJD(1);
t19 = sin(qJ(2));
t27 = t19 * t32;
t7 = t14 * qJ(3) + t27;
t31 = t7 ^ 2 / 0.2e1;
t30 = qJD(4) * t18;
t29 = qJD(4) * t20;
t13 = t14 ^ 2;
t28 = t20 * t13 * t18;
t21 = cos(qJ(2));
t26 = t21 * t32;
t25 = t14 * t30;
t24 = qJD(3) - t26;
t22 = qJD(1) ^ 2;
t15 = qJD(4) ^ 2 / 0.2e1;
t12 = t13 / 0.2e1;
t11 = t14 * t29;
t10 = t13 * t35;
t9 = t13 * t36;
t6 = -t14 * pkin(2) + t24;
t5 = (-pkin(2) - pkin(7)) * t14 + t24;
t3 = qJD(4) * qJ(5) + t18 * t5;
t2 = -qJD(4) * pkin(4) - t20 * t5 + qJD(5);
t1 = t27 + (pkin(4) * t18 - qJ(5) * t20 + qJ(3)) * t14;
t4 = [0, 0, 0, 0, 0, t22 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t14 * t26, -t14 * t27, 0, (t19 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t22, t12, 0, 0, 0, 0, 0, 0, t6 * t14, t7 * t14, t31 + t6 ^ 2 / 0.2e1, t10, -t28, t11, t9, -t25, t15, t5 * t29 + t7 * t34, -t5 * t30 + t7 * t33, (-t16 - t17) * t5 * t14, t31 + (t36 + t35) * t5 ^ 2, t10, t11, t28, t15, t25, t9, -t2 * qJD(4) + t1 * t34, (-t18 * t3 + t2 * t20) * t14, t3 * qJD(4) - t1 * t33, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t4;
