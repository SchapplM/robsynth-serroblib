% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:37
% EndTime: 2019-12-31 17:56:38
% DurationCPUTime: 0.10s
% Computational Cost: add. (121->26), mult. (227->67), div. (0->0), fcn. (82->6), ass. (0->26)
t13 = -qJD(1) + qJD(4);
t12 = t13 ^ 2;
t28 = t12 / 0.2e1;
t16 = sin(pkin(8));
t11 = (pkin(1) * t16 + qJ(3)) * qJD(1);
t19 = sin(qJ(4));
t21 = cos(qJ(4));
t17 = cos(pkin(8));
t24 = -pkin(1) * t17 - pkin(2);
t8 = qJD(3) + (-pkin(3) + t24) * qJD(1);
t6 = t21 * t11 + t19 * t8;
t22 = qJD(1) ^ 2;
t27 = pkin(1) * t22;
t5 = -t19 * t11 + t21 * t8;
t3 = -t13 * pkin(4) - t5;
t26 = t13 * t3;
t25 = qJD(5) * t13;
t20 = cos(qJ(5));
t18 = sin(qJ(5));
t15 = t22 / 0.2e1;
t14 = qJD(2) ^ 2 / 0.2e1;
t10 = t24 * qJD(1) + qJD(3);
t4 = t13 * pkin(7) + t6;
t2 = -t18 * qJD(2) + t20 * t4;
t1 = -t20 * qJD(2) - t18 * t4;
t7 = [0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t17 * t27, -t16 * t27, 0, t14 + (t16 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t22, 0, 0, 0, t15, 0, 0, -t10 * qJD(1), 0, t11 * qJD(1), t11 ^ 2 / 0.2e1 + t14 + t10 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t28, t5 * t13, -t6 * t13, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14, t18 ^ 2 * t28, t18 * t12 * t20, t18 * t25, t20 ^ 2 * t28, t20 * t25, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t20 * t26, -t2 * qJD(5) + t18 * t26, (-t1 * t18 + t2 * t20) * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
