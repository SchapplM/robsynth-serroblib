% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:53
% EndTime: 2019-12-31 16:38:53
% DurationCPUTime: 0.08s
% Computational Cost: add. (46->20), mult. (124->49), div. (0->0), fcn. (39->4), ass. (0->18)
t13 = qJD(1) ^ 2;
t8 = t13 / 0.2e1;
t19 = pkin(1) * t13;
t9 = sin(pkin(6));
t5 = (-pkin(1) * t9 - qJ(3)) * qJD(1);
t18 = t5 * qJD(1);
t17 = t5 ^ 2 / 0.2e1;
t16 = qJD(1) * qJD(4);
t10 = cos(pkin(6));
t15 = -pkin(1) * t10 - pkin(2);
t12 = cos(qJ(4));
t11 = sin(qJ(4));
t7 = qJD(2) ^ 2 / 0.2e1;
t4 = t15 * qJD(1) + qJD(3);
t3 = qJD(3) + (-pkin(5) + t15) * qJD(1);
t2 = qJD(2) * t12 + t11 * t3;
t1 = -qJD(2) * t11 + t12 * t3;
t6 = [0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10 * t19, -t9 * t19, 0, t7 + (t9 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t13, t8, 0, 0, 0, 0, 0, 0, t4 * qJD(1), -t18, t7 + t17 + t4 ^ 2 / 0.2e1, t12 ^ 2 * t8, -t12 * t13 * t11, t12 * t16, t11 ^ 2 * t8, -t11 * t16, qJD(4) ^ 2 / 0.2e1, qJD(4) * t1 - t11 * t18, -qJD(4) * t2 - t12 * t18, (-t1 * t12 - t11 * t2) * qJD(1), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t17;];
T_reg = t6;
