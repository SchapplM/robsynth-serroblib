% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:11
% EndTime: 2019-12-31 17:41:11
% DurationCPUTime: 0.10s
% Computational Cost: add. (82->32), mult. (226->68), div. (0->0), fcn. (87->2), ass. (0->28)
t21 = qJD(2) ^ 2;
t31 = t21 / 0.2e1;
t30 = pkin(3) + pkin(4);
t20 = cos(qJ(3));
t29 = t20 * t21;
t19 = sin(qJ(3));
t27 = qJD(2) * t20;
t8 = pkin(6) * t27 + t19 * qJD(1);
t28 = qJD(2) * t19;
t26 = qJ(5) * qJD(2);
t25 = qJD(2) * qJD(3);
t24 = t19 * t29;
t6 = qJD(3) * qJ(4) + t8;
t23 = qJ(4) * t19 + pkin(2);
t7 = -pkin(6) * t28 + t20 * qJD(1);
t22 = qJD(4) - t7;
t17 = qJD(1) ^ 2 / 0.2e1;
t16 = qJD(3) ^ 2 / 0.2e1;
t12 = t20 * t25;
t11 = t19 * t25;
t10 = t20 ^ 2 * t31;
t9 = t19 ^ 2 * t31;
t5 = (-pkin(3) * t20 - t23) * qJD(2);
t4 = -qJD(3) * pkin(3) + t22;
t3 = -t20 * t26 + t6;
t2 = qJD(5) + (t30 * t20 + t23) * qJD(2);
t1 = -t30 * qJD(3) - t19 * t26 + t22;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, t31, 0, 0, 0, t17, t9, t24, t11, t10, t12, t16, pkin(2) * t29 + t7 * qJD(3), -t21 * pkin(2) * t19 - t8 * qJD(3), (-t19 * t7 + t20 * t8) * qJD(2), t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + pkin(2) ^ 2 * t31, t9, t11, -t24, t16, -t12, t10, -t4 * qJD(3) - t5 * t27, (t19 * t4 + t20 * t6) * qJD(2), t6 * qJD(3) - t5 * t28, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t9, -t24, -t11, t10, t12, t16, -t1 * qJD(3) + t2 * t27, t3 * qJD(3) + t2 * t28, (-t1 * t19 - t20 * t3) * qJD(2), t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t13;
