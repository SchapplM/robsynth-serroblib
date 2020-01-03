% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:17
% EndTime: 2019-12-31 16:59:17
% DurationCPUTime: 0.10s
% Computational Cost: add. (64->30), mult. (208->67), div. (0->0), fcn. (67->2), ass. (0->27)
t19 = qJD(1) ^ 2;
t30 = t19 / 0.2e1;
t29 = pkin(2) + pkin(3);
t18 = cos(qJ(2));
t28 = t18 * t19;
t26 = qJD(1) * t18;
t6 = pkin(5) * t26 + qJD(2) * qJ(3);
t17 = sin(qJ(2));
t27 = qJD(1) * t17;
t25 = pkin(5) * t27 + qJD(3);
t24 = qJ(4) * qJD(1);
t23 = qJD(1) * qJD(2);
t22 = t17 * t28;
t9 = t17 * t23;
t10 = t18 * t23;
t21 = qJ(3) * t17 + pkin(1);
t16 = t18 ^ 2;
t15 = t17 ^ 2;
t13 = qJD(2) ^ 2 / 0.2e1;
t8 = t16 * t30;
t7 = t15 * t30;
t5 = -qJD(2) * pkin(2) + t25;
t4 = (-pkin(2) * t18 - t21) * qJD(1);
t3 = -t18 * t24 + t6;
t2 = -t29 * qJD(2) - t17 * t24 + t25;
t1 = qJD(4) + (t29 * t18 + t21) * qJD(1);
t11 = [0, 0, 0, 0, 0, t30, 0, 0, 0, 0, t7, t22, t9, t8, t10, t13, pkin(1) * t28 - pkin(5) * t9, -t19 * pkin(1) * t17 - pkin(5) * t10, (t15 + t16) * t19 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t16 / 0.2e1 + t15 / 0.2e1) * pkin(5) ^ 2) * t19, t7, t9, -t22, t13, -t10, t8, -t5 * qJD(2) - t4 * t26, (t17 * t5 + t18 * t6) * qJD(1), t6 * qJD(2) - t4 * t27, t6 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t7, -t22, -t9, t8, t10, t13, -t2 * qJD(2) + t1 * t26, t3 * qJD(2) + t1 * t27, (-t17 * t2 - t18 * t3) * qJD(1), t3 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t11;
