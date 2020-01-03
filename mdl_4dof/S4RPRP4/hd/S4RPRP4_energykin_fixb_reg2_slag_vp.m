% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:52
% EndTime: 2019-12-31 16:43:52
% DurationCPUTime: 0.10s
% Computational Cost: add. (60->24), mult. (186->62), div. (0->0), fcn. (72->4), ass. (0->25)
t19 = qJD(1) ^ 2;
t14 = t19 / 0.2e1;
t27 = pkin(1) * t19;
t17 = sin(qJ(3));
t18 = cos(qJ(3));
t15 = sin(pkin(6));
t7 = (pkin(1) * t15 + pkin(5)) * qJD(1);
t5 = t17 * qJD(2) + t18 * t7;
t26 = qJD(1) * t17;
t25 = qJD(1) * t18;
t24 = qJD(1) * qJD(3);
t23 = t17 * t19 * t18;
t16 = cos(pkin(6));
t22 = -pkin(1) * t16 - pkin(2);
t21 = t18 * t24;
t4 = t18 * qJD(2) - t17 * t7;
t13 = qJD(3) ^ 2 / 0.2e1;
t11 = t17 * t24;
t10 = t18 ^ 2 * t14;
t9 = t17 ^ 2 * t14;
t8 = t22 * qJD(1);
t3 = (-pkin(3) * t18 - qJ(4) * t17 + t22) * qJD(1);
t2 = qJD(3) * qJ(4) + t5;
t1 = -qJD(3) * pkin(3) + qJD(4) - t4;
t6 = [0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16 * t27, -t15 * t27, 0, qJD(2) ^ 2 / 0.2e1 + (t15 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t19, t9, t23, t11, t10, t21, t13, t4 * qJD(3) - t8 * t25, -t5 * qJD(3) + t8 * t26, (-t17 * t4 + t18 * t5) * qJD(1), t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t9, t11, -t23, t13, -t21, t10, -t1 * qJD(3) - t3 * t25, (t1 * t17 + t18 * t2) * qJD(1), t2 * qJD(3) - t3 * t26, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
