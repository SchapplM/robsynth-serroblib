% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:53
% EndTime: 2019-12-31 18:05:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (113->31), mult. (228->77), div. (0->0), fcn. (82->4), ass. (0->28)
t24 = qJD(1) ^ 2;
t17 = t24 / 0.2e1;
t32 = cos(qJ(5));
t31 = -pkin(1) - qJ(3);
t23 = cos(qJ(4));
t30 = qJD(1) * t23;
t14 = qJD(1) * qJ(2) + qJD(3);
t10 = -qJD(1) * pkin(6) + t14;
t29 = qJD(4) * t10;
t11 = -t31 * qJD(1) - qJD(2);
t28 = t11 * qJD(1);
t22 = sin(qJ(4));
t27 = t22 * qJD(1);
t26 = t11 ^ 2 / 0.2e1;
t25 = qJD(1) * qJD(4);
t21 = sin(qJ(5));
t20 = t23 ^ 2;
t19 = t22 ^ 2;
t15 = -qJD(1) * pkin(1) + qJD(2);
t13 = qJD(5) + t27;
t8 = t21 * qJD(4) + t32 * t30;
t6 = -t32 * qJD(4) + t21 * t30;
t5 = -qJD(4) * pkin(4) - t23 * t10;
t4 = qJD(4) * pkin(7) + t22 * t10;
t3 = -qJD(2) + (pkin(4) * t22 - pkin(7) * t23 - t31) * qJD(1);
t2 = t21 * t3 + t32 * t4;
t1 = -t21 * t4 + t32 * t3;
t7 = [0, 0, 0, 0, 0, t17, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, t15 * qJD(1), t24 * qJ(2), qJ(2) ^ 2 * t17 + t15 ^ 2 / 0.2e1, t17, 0, 0, 0, 0, 0, 0, t14 * qJD(1), t28, t26 + t14 ^ 2 / 0.2e1, t20 * t17, -t23 * t24 * t22, t23 * t25, t19 * t17, -t22 * t25, qJD(4) ^ 2 / 0.2e1, t11 * t27 + t23 * t29, -t22 * t29 + t23 * t28, (-t19 - t20) * t10 * qJD(1), t26 + (t19 / 0.2e1 + t20 / 0.2e1) * t10 ^ 2, t8 ^ 2 / 0.2e1, -t8 * t6, t8 * t13, t6 ^ 2 / 0.2e1, -t6 * t13, t13 ^ 2 / 0.2e1, t1 * t13 + t5 * t6, -t2 * t13 + t5 * t8, -t1 * t8 - t2 * t6, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
