% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:23
% EndTime: 2019-12-31 17:32:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (88->29), mult. (226->75), div. (0->0), fcn. (133->6), ass. (0->27)
t21 = qJD(3) ^ 2;
t29 = t21 / 0.2e1;
t28 = cos(qJ(5));
t16 = sin(pkin(8));
t27 = qJD(3) * t16;
t17 = cos(pkin(8));
t26 = qJD(3) * t17;
t25 = t17 * qJD(1);
t24 = qJD(2) * qJD(3);
t19 = sin(qJ(3));
t13 = qJD(3) * qJ(4) + t19 * qJD(2);
t6 = -t16 * qJD(1) + t17 * t13;
t20 = cos(qJ(3));
t23 = -t20 * qJD(2) + qJD(4);
t22 = qJD(2) ^ 2;
t18 = sin(qJ(5));
t15 = qJD(1) ^ 2 / 0.2e1;
t12 = -qJD(3) * pkin(3) + t23;
t10 = (-pkin(4) * t17 - pkin(3)) * qJD(3) + t23;
t9 = (t28 * t16 + t17 * t18) * qJD(3);
t7 = t18 * t27 - t28 * t26;
t5 = -t16 * t13 - t25;
t4 = pkin(6) * t26 + t6;
t3 = -t25 + (-pkin(6) * qJD(3) - t13) * t16;
t2 = t18 * t3 + t28 * t4;
t1 = -t18 * t4 + t28 * t3;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 + t22 / 0.2e1, 0, 0, 0, 0, 0, t29, t20 * t24, -t19 * t24, 0, t15 + (t19 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1) * t22, t16 ^ 2 * t29, t16 * t21 * t17, 0, t17 ^ 2 * t29, 0, 0, -t12 * t26, t12 * t27, (-t16 * t5 + t17 * t6) * qJD(3), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * qJD(5), t7 ^ 2 / 0.2e1, -t7 * qJD(5), qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t10 * t7, -t2 * qJD(5) + t10 * t9, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
