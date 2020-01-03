% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:13
% EndTime: 2019-12-31 17:22:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (102->16), mult. (171->56), div. (0->0), fcn. (70->6), ass. (0->26)
t14 = sin(qJ(4));
t12 = t14 ^ 2;
t29 = t12 / 0.2e1;
t17 = cos(qJ(4));
t13 = t17 ^ 2;
t28 = t13 / 0.2e1;
t15 = sin(qJ(3));
t18 = cos(qJ(3));
t16 = sin(qJ(2));
t26 = pkin(1) * qJD(1);
t23 = t16 * t26;
t11 = qJD(1) + qJD(2);
t19 = cos(qJ(2));
t22 = t19 * t26;
t7 = pkin(2) * t11 + t22;
t5 = t15 * t7 + t18 * t23;
t10 = qJD(3) + t11;
t4 = -t15 * t23 + t18 * t7;
t2 = -pkin(3) * t10 - t4;
t27 = t10 * t2;
t25 = qJD(4) * t14;
t24 = qJD(4) * t17;
t20 = qJD(1) ^ 2;
t9 = t10 ^ 2;
t3 = pkin(7) * t10 + t5;
t1 = [0, 0, 0, 0, 0, t20 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2 / 0.2e1, t11 * t22, -t11 * t23, 0, (t16 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t20, 0, 0, 0, 0, 0, t9 / 0.2e1, t4 * t10, -t5 * t10, 0, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t9 * t29, t14 * t9 * t17, t10 * t25, t9 * t28, t10 * t24, qJD(4) ^ 2 / 0.2e1, -t17 * t27 - t3 * t25, t14 * t27 - t3 * t24, (t12 + t13) * t3 * t10, t2 ^ 2 / 0.2e1 + (t28 + t29) * t3 ^ 2;];
T_reg = t1;
