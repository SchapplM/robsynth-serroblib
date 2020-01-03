% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:54
% EndTime: 2019-12-31 17:11:55
% DurationCPUTime: 0.11s
% Computational Cost: add. (101->35), mult. (282->82), div. (0->0), fcn. (119->4), ass. (0->34)
t24 = qJD(1) ^ 2;
t36 = t24 / 0.2e1;
t35 = -pkin(2) - pkin(6);
t23 = cos(qJ(2));
t34 = t23 * t24;
t33 = qJD(1) * t23;
t21 = sin(qJ(2));
t32 = t21 * qJD(1);
t31 = pkin(5) * t32 + qJD(3);
t30 = qJD(2) * qJ(3);
t29 = qJD(1) * qJD(2);
t28 = t21 * t29;
t27 = t23 * t29;
t26 = -qJ(3) * t21 - pkin(1);
t22 = cos(qJ(4));
t20 = sin(qJ(4));
t19 = t23 ^ 2;
t18 = t21 ^ 2;
t17 = qJD(2) ^ 2 / 0.2e1;
t15 = t19 * t36;
t14 = t18 * t36;
t13 = qJD(4) + t32;
t12 = t21 * t34;
t11 = -pkin(5) * t33 - t30;
t10 = -qJD(2) * pkin(2) + t31;
t9 = t22 * qJD(2) - t20 * t33;
t7 = t20 * qJD(2) + t22 * t33;
t6 = (-pkin(2) * t23 + t26) * qJD(1);
t5 = t30 + (pkin(3) + pkin(5)) * t33;
t4 = pkin(3) * t32 + t35 * qJD(2) + t31;
t3 = (t35 * t23 + t26) * qJD(1);
t2 = t20 * t4 + t22 * t3;
t1 = -t20 * t3 + t22 * t4;
t8 = [0, 0, 0, 0, 0, t36, 0, 0, 0, 0, t14, t12, t28, t15, t27, t17, pkin(1) * t34 - pkin(5) * t28, -t24 * pkin(1) * t21 - pkin(5) * t27, (t18 + t19) * t24 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t19 / 0.2e1 + t18 / 0.2e1) * pkin(5) ^ 2) * t24, t17, -t28, -t27, t14, t12, t15, (t10 * t21 - t11 * t23) * qJD(1), t10 * qJD(2) + t6 * t33, -t11 * qJD(2) - t6 * t32, t6 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t13, t7 ^ 2 / 0.2e1, -t7 * t13, t13 ^ 2 / 0.2e1, t1 * t13 + t5 * t7, -t2 * t13 + t5 * t9, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t8;
