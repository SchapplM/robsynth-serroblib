% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:50
% EndTime: 2019-12-31 16:57:50
% DurationCPUTime: 0.11s
% Computational Cost: add. (116->30), mult. (363->75), div. (0->0), fcn. (196->4), ass. (0->32)
t23 = qJD(1) ^ 2;
t36 = t23 / 0.2e1;
t20 = sin(pkin(6));
t21 = sin(qJ(2));
t22 = cos(qJ(2));
t30 = cos(pkin(6));
t11 = (t20 * t22 + t30 * t21) * qJD(1);
t28 = qJD(1) * t22;
t29 = qJD(1) * t21;
t9 = t20 * t29 - t30 * t28;
t35 = t11 * t9;
t33 = pkin(5) + qJ(3);
t13 = qJD(2) * pkin(2) - t33 * t29;
t14 = t33 * t28;
t5 = t20 * t13 + t30 * t14;
t34 = t22 * t23;
t32 = qJD(2) * t9;
t31 = t9 ^ 2 / 0.2e1;
t27 = qJD(1) * qJD(2);
t26 = t21 * t27;
t25 = t22 * t27;
t4 = t30 * t13 - t20 * t14;
t15 = qJD(3) + (-pkin(2) * t22 - pkin(1)) * qJD(1);
t19 = t22 ^ 2;
t18 = t21 ^ 2;
t17 = qJD(2) ^ 2 / 0.2e1;
t8 = t11 * qJD(2);
t7 = t11 ^ 2 / 0.2e1;
t3 = qJD(2) * qJ(4) + t5;
t2 = -qJD(2) * pkin(3) + qJD(4) - t4;
t1 = t9 * pkin(3) - t11 * qJ(4) + t15;
t6 = [0, 0, 0, 0, 0, t36, 0, 0, 0, 0, t18 * t36, t21 * t34, t26, t19 * t36, t25, t17, pkin(1) * t34 - pkin(5) * t26, -t23 * pkin(1) * t21 - pkin(5) * t25, (t18 + t19) * t23 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t19 / 0.2e1 + t18 / 0.2e1) * pkin(5) ^ 2) * t23, t7, -t35, t8, t31, -t32, t17, t4 * qJD(2) + t15 * t9, -t5 * qJD(2) + t15 * t11, -t4 * t11 - t5 * t9, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t7, t8, t35, t17, t32, t31, -t2 * qJD(2) + t1 * t9, t2 * t11 - t3 * t9, t3 * qJD(2) - t1 * t11, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t6;
