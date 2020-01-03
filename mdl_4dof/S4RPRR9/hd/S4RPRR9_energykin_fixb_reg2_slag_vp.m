% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRR9
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:28
% EndTime: 2019-12-31 16:56:28
% DurationCPUTime: 0.09s
% Computational Cost: add. (87->27), mult. (206->69), div. (0->0), fcn. (82->4), ass. (0->24)
t21 = qJD(1) ^ 2;
t15 = t21 / 0.2e1;
t26 = cos(qJ(4));
t25 = t21 * qJ(2);
t20 = cos(qJ(3));
t24 = qJD(1) * t20;
t10 = qJD(2) + (-pkin(1) - pkin(5)) * qJD(1);
t23 = qJD(3) * t10;
t22 = qJD(1) * qJD(3);
t19 = sin(qJ(3));
t18 = sin(qJ(4));
t17 = t20 ^ 2;
t16 = t19 ^ 2;
t13 = qJ(2) ^ 2 * t15;
t12 = -qJD(1) * pkin(1) + qJD(2);
t11 = t19 * qJD(1) + qJD(4);
t8 = t18 * qJD(3) + t26 * t24;
t6 = -t26 * qJD(3) + t18 * t24;
t5 = -qJD(3) * pkin(3) - t20 * t10;
t4 = qJD(3) * pkin(6) + t19 * t10;
t3 = (pkin(3) * t19 - pkin(6) * t20 + qJ(2)) * qJD(1);
t2 = t18 * t3 + t26 * t4;
t1 = -t18 * t4 + t26 * t3;
t7 = [0, 0, 0, 0, 0, t15, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, t12 * qJD(1), t25, t13 + t12 ^ 2 / 0.2e1, t17 * t15, -t20 * t21 * t19, t20 * t22, t16 * t15, -t19 * t22, qJD(3) ^ 2 / 0.2e1, t19 * t25 + t20 * t23, -t19 * t23 + t20 * t25, (-t16 - t17) * t10 * qJD(1), t13 + (t16 / 0.2e1 + t17 / 0.2e1) * t10 ^ 2, t8 ^ 2 / 0.2e1, -t8 * t6, t8 * t11, t6 ^ 2 / 0.2e1, -t6 * t11, t11 ^ 2 / 0.2e1, t1 * t11 + t5 * t6, -t2 * t11 + t5 * t8, -t1 * t8 - t2 * t6, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
