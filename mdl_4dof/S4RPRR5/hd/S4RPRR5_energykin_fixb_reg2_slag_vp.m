% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:39
% EndTime: 2019-12-31 16:51:39
% DurationCPUTime: 0.10s
% Computational Cost: add. (80->17), mult. (138->50), div. (0->0), fcn. (40->4), ass. (0->23)
t15 = sin(qJ(4));
t13 = t15 ^ 2;
t25 = t13 / 0.2e1;
t17 = cos(qJ(4));
t14 = t17 ^ 2;
t24 = t14 / 0.2e1;
t19 = qJD(1) ^ 2;
t12 = t19 / 0.2e1;
t16 = sin(qJ(3));
t18 = cos(qJ(3));
t20 = qJ(2) * qJD(1);
t7 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t5 = t16 * t7 + t18 * t20;
t11 = -qJD(1) + qJD(3);
t4 = -t16 * t20 + t18 * t7;
t2 = -pkin(3) * t11 - t4;
t23 = t11 * t2;
t22 = qJD(4) * t15;
t21 = qJD(4) * t17;
t10 = t11 ^ 2;
t9 = -qJD(1) * pkin(1) + qJD(2);
t3 = pkin(6) * t11 + t5;
t1 = [0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, -t9 * qJD(1), 0, t19 * qJ(2), qJ(2) ^ 2 * t12 + t9 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t10 / 0.2e1, t4 * t11, -t5 * t11, 0, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t10 * t25, t15 * t10 * t17, t11 * t22, t10 * t24, t11 * t21, qJD(4) ^ 2 / 0.2e1, -t17 * t23 - t3 * t22, t15 * t23 - t3 * t21, (t13 + t14) * t3 * t11, t2 ^ 2 / 0.2e1 + (t24 + t25) * t3 ^ 2;];
T_reg = t1;
