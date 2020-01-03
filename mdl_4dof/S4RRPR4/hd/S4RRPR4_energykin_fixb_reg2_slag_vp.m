% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:32
% EndTime: 2019-12-31 17:02:32
% DurationCPUTime: 0.10s
% Computational Cost: add. (137->24), mult. (234->71), div. (0->0), fcn. (115->6), ass. (0->30)
t17 = sin(pkin(7));
t14 = t17 ^ 2;
t33 = t14 / 0.2e1;
t18 = cos(pkin(7));
t15 = t18 ^ 2;
t32 = t15 / 0.2e1;
t31 = cos(qJ(4));
t16 = qJD(1) + qJD(2);
t30 = t16 * t17;
t29 = t16 * t18;
t28 = pkin(1) * qJD(1);
t20 = sin(qJ(2));
t27 = t20 * t28;
t21 = cos(qJ(2));
t26 = t21 * t28;
t11 = t16 * qJ(3) + t27;
t25 = pkin(6) * t16 + t11;
t24 = qJD(3) - t26;
t22 = qJD(1) ^ 2;
t19 = sin(qJ(4));
t13 = t16 ^ 2;
t9 = -t16 * pkin(2) + t24;
t8 = (t31 * t17 + t18 * t19) * t16;
t6 = t19 * t30 - t31 * t29;
t5 = (-pkin(3) * t18 - pkin(2)) * t16 + t24;
t4 = t25 * t18;
t3 = t25 * t17;
t2 = -t19 * t3 + t31 * t4;
t1 = -t19 * t4 - t31 * t3;
t7 = [0, 0, 0, 0, 0, t22 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 / 0.2e1, t16 * t26, -t16 * t27, 0, (t20 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t22, t13 * t33, t17 * t13 * t18, 0, t13 * t32, 0, 0, -t9 * t29, t9 * t30, (t14 + t15) * t16 * t11, t9 ^ 2 / 0.2e1 + (t32 + t33) * t11 ^ 2, t8 ^ 2 / 0.2e1, -t8 * t6, t8 * qJD(4), t6 ^ 2 / 0.2e1, -t6 * qJD(4), qJD(4) ^ 2 / 0.2e1, t1 * qJD(4) + t5 * t6, -t2 * qJD(4) + t5 * t8, -t1 * t8 - t2 * t6, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
