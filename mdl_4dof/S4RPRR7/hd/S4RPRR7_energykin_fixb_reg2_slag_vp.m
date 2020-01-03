% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:04
% EndTime: 2019-12-31 16:54:04
% DurationCPUTime: 0.12s
% Computational Cost: add. (159->35), mult. (462->89), div. (0->0), fcn. (289->6), ass. (0->30)
t29 = qJD(1) ^ 2;
t36 = t29 / 0.2e1;
t35 = cos(qJ(3));
t34 = cos(qJ(4));
t33 = pkin(5) + qJ(2);
t25 = sin(pkin(7));
t32 = qJD(1) * t25;
t17 = t33 * t32;
t26 = cos(pkin(7));
t31 = qJD(1) * t26;
t18 = t33 * t31;
t28 = sin(qJ(3));
t7 = -t28 * t17 + t35 * t18;
t14 = t28 * t32 - t35 * t31;
t6 = -t35 * t17 - t28 * t18;
t19 = qJD(2) + (-pkin(2) * t26 - pkin(1)) * qJD(1);
t27 = sin(qJ(4));
t24 = t26 ^ 2;
t23 = t25 ^ 2;
t21 = -qJD(1) * pkin(1) + qJD(2);
t16 = (t35 * t25 + t26 * t28) * qJD(1);
t11 = qJD(4) + t14;
t10 = t27 * qJD(3) + t34 * t16;
t8 = -t34 * qJD(3) + t27 * t16;
t5 = qJD(3) * pkin(6) + t7;
t4 = -qJD(3) * pkin(3) - t6;
t3 = t14 * pkin(3) - t16 * pkin(6) + t19;
t2 = t27 * t3 + t34 * t5;
t1 = -t27 * t5 + t34 * t3;
t9 = [0, 0, 0, 0, 0, t36, 0, 0, 0, 0, t23 * t36, t25 * t29 * t26, 0, t24 * t36, 0, 0, -t21 * t31, t21 * t32, (t23 + t24) * t29 * qJ(2), t21 ^ 2 / 0.2e1 + (t24 / 0.2e1 + t23 / 0.2e1) * qJ(2) ^ 2 * t29, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * qJD(3), t14 ^ 2 / 0.2e1, -t14 * qJD(3), qJD(3) ^ 2 / 0.2e1, t6 * qJD(3) + t19 * t14, -t7 * qJD(3) + t19 * t16, -t7 * t14 - t6 * t16, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t11, t8 ^ 2 / 0.2e1, -t8 * t11, t11 ^ 2 / 0.2e1, t1 * t11 + t4 * t8, t4 * t10 - t2 * t11, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg = t9;
