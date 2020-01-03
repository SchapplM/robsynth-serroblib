% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:08
% EndTime: 2019-12-31 17:13:09
% DurationCPUTime: 0.09s
% Computational Cost: add. (93->20), mult. (177->61), div. (0->0), fcn. (61->4), ass. (0->31)
t17 = sin(qJ(3));
t15 = t17 ^ 2;
t32 = t15 / 0.2e1;
t19 = cos(qJ(3));
t16 = t19 ^ 2;
t31 = t16 / 0.2e1;
t13 = qJD(1) + qJD(2);
t30 = t13 * t17;
t29 = t13 * t19;
t28 = pkin(1) * qJD(1);
t27 = qJD(3) * t17;
t26 = qJD(3) * t19;
t18 = sin(qJ(2));
t25 = t18 * t28;
t20 = cos(qJ(2));
t24 = t20 * t28;
t5 = t13 * pkin(6) + t25;
t23 = qJ(4) * t13 + t5;
t21 = qJD(1) ^ 2;
t14 = qJD(3) ^ 2 / 0.2e1;
t12 = t13 ^ 2;
t11 = t13 * t26;
t10 = t13 * t27;
t9 = t12 * t31;
t8 = t12 * t32;
t7 = t17 * t12 * t19;
t6 = -t13 * pkin(2) - t24;
t3 = -t24 + qJD(4) + (-pkin(3) * t19 - pkin(2)) * t13;
t2 = t23 * t19;
t1 = qJD(3) * pkin(3) - t23 * t17;
t4 = [0, 0, 0, 0, 0, t21 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 / 0.2e1, t13 * t24, -t13 * t25, 0, (t18 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t21, t8, t7, t10, t9, t11, t14, -t5 * t27 - t6 * t29, -t5 * t26 + t6 * t30, (t15 + t16) * t5 * t13, t6 ^ 2 / 0.2e1 + (t31 + t32) * t5 ^ 2, t8, t7, t10, t9, t11, t14, t1 * qJD(3) - t3 * t29, -t2 * qJD(3) + t3 * t30, (-t1 * t17 + t19 * t2) * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t4;
