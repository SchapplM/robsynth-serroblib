% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:38
% EndTime: 2019-12-31 16:32:39
% DurationCPUTime: 0.09s
% Computational Cost: add. (70->24), mult. (209->69), div. (0->0), fcn. (113->4), ass. (0->23)
t20 = qJD(2) ^ 2;
t26 = t20 / 0.2e1;
t25 = cos(qJ(4));
t19 = cos(qJ(3));
t24 = t19 * t20;
t18 = sin(qJ(3));
t22 = qJD(2) * t19;
t9 = pkin(5) * t22 + t18 * qJD(1);
t23 = qJD(2) * t18;
t21 = qJD(2) * qJD(3);
t17 = sin(qJ(4));
t16 = qJD(1) ^ 2 / 0.2e1;
t15 = qJD(3) + qJD(4);
t14 = t19 * qJD(1);
t10 = (-pkin(3) * t19 - pkin(2)) * qJD(2);
t8 = -pkin(5) * t23 + t14;
t7 = (t17 * t19 + t25 * t18) * qJD(2);
t5 = t17 * t23 - t25 * t22;
t4 = pkin(6) * t22 + t9;
t3 = qJD(3) * pkin(3) + t14 + (-pkin(6) - pkin(5)) * t23;
t2 = t17 * t3 + t25 * t4;
t1 = -t17 * t4 + t25 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, t26, 0, 0, 0, t16, t18 ^ 2 * t26, t18 * t24, t18 * t21, t19 ^ 2 * t26, t19 * t21, qJD(3) ^ 2 / 0.2e1, pkin(2) * t24 + t8 * qJD(3), -t20 * pkin(2) * t18 - t9 * qJD(3), (-t18 * t8 + t19 * t9) * qJD(2), t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + pkin(2) ^ 2 * t26, t7 ^ 2 / 0.2e1, -t7 * t5, t7 * t15, t5 ^ 2 / 0.2e1, -t5 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t10 * t5, t10 * t7 - t2 * t15, -t1 * t7 - t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t6;
