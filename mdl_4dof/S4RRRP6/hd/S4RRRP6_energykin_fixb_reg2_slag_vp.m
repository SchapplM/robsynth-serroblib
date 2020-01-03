% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:09
% EndTime: 2019-12-31 17:19:09
% DurationCPUTime: 0.11s
% Computational Cost: add. (126->33), mult. (347->78), div. (0->0), fcn. (180->4), ass. (0->32)
t30 = qJD(1) ^ 2;
t39 = t30 / 0.2e1;
t28 = sin(qJ(2));
t29 = cos(qJ(2));
t11 = (-pkin(2) * t29 - pkin(6) * t28 - pkin(1)) * qJD(1);
t35 = t29 * qJD(1);
t19 = pkin(5) * t35 + qJD(2) * pkin(6);
t27 = sin(qJ(3));
t38 = cos(qJ(3));
t4 = t27 * t11 + t38 * t19;
t37 = t29 * t30;
t36 = qJD(1) * t28;
t34 = qJD(1) * qJD(2);
t33 = t28 * t34;
t32 = t29 * t34;
t3 = t38 * t11 - t27 * t19;
t18 = -qJD(2) * pkin(2) + pkin(5) * t36;
t26 = t29 ^ 2;
t25 = t28 ^ 2;
t21 = -qJD(3) + t35;
t20 = t21 ^ 2 / 0.2e1;
t16 = t27 * qJD(2) + t38 * t36;
t14 = -t38 * qJD(2) + t27 * t36;
t13 = t16 ^ 2 / 0.2e1;
t12 = t14 ^ 2 / 0.2e1;
t8 = t16 * t21;
t7 = t14 * t21;
t6 = t14 * pkin(3) + qJD(4) + t18;
t5 = t16 * t14;
t2 = -t14 * qJ(4) + t4;
t1 = -t21 * pkin(3) - t16 * qJ(4) + t3;
t9 = [0, 0, 0, 0, 0, t39, 0, 0, 0, 0, t25 * t39, t28 * t37, t33, t26 * t39, t32, qJD(2) ^ 2 / 0.2e1, pkin(1) * t37 - pkin(5) * t33, -t30 * pkin(1) * t28 - pkin(5) * t32, (t25 + t26) * t30 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t26 / 0.2e1 + t25 / 0.2e1) * pkin(5) ^ 2) * t30, t13, -t5, -t8, t12, t7, t20, t18 * t14 - t3 * t21, t18 * t16 + t4 * t21, -t4 * t14 - t3 * t16, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t13, -t5, -t8, t12, t7, t20, -t1 * t21 + t6 * t14, t6 * t16 + t2 * t21, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg = t9;
