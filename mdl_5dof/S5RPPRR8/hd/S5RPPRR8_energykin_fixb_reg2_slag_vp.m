% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:13
% EndTime: 2019-12-31 18:01:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (167->27), mult. (278->67), div. (0->0), fcn. (108->6), ass. (0->28)
t16 = -qJD(1) + qJD(4);
t15 = t16 ^ 2;
t29 = t15 / 0.2e1;
t25 = qJD(1) ^ 2;
t18 = t25 / 0.2e1;
t13 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t26 = qJ(2) * qJD(1);
t11 = t19 * t13 + t20 * t26;
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t12 = t20 * t13;
t9 = t12 + (-qJ(2) * t19 - pkin(3)) * qJD(1);
t6 = t24 * t11 + t22 * t9;
t5 = -t22 * t11 + t24 * t9;
t3 = -t16 * pkin(4) - t5;
t28 = t16 * t3;
t27 = qJD(5) * t16;
t23 = cos(qJ(5));
t21 = sin(qJ(5));
t17 = qJD(3) ^ 2 / 0.2e1;
t14 = -pkin(1) * qJD(1) + qJD(2);
t10 = -t19 * t26 + t12;
t4 = t16 * pkin(7) + t6;
t2 = t21 * qJD(3) + t23 * t4;
t1 = t23 * qJD(3) - t21 * t4;
t7 = [0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, -t14 * qJD(1), 0, t25 * qJ(2), qJ(2) ^ 2 * t18 + t14 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t18, -t10 * qJD(1), t11 * qJD(1), 0, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t17, 0, 0, 0, 0, 0, t29, t5 * t16, -t6 * t16, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17, t21 ^ 2 * t29, t21 * t15 * t23, t21 * t27, t23 ^ 2 * t29, t23 * t27, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t23 * t28, -t2 * qJD(5) + t21 * t28, (-t1 * t21 + t2 * t23) * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
