% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:06
% EndTime: 2019-12-31 17:58:06
% DurationCPUTime: 0.16s
% Computational Cost: add. (242->44), mult. (606->109), div. (0->0), fcn. (375->8), ass. (0->35)
t36 = qJD(1) ^ 2;
t29 = t36 / 0.2e1;
t31 = sin(pkin(8));
t24 = (pkin(1) * t31 + qJ(3)) * qJD(1);
t32 = cos(pkin(9));
t27 = t32 * qJD(2);
t30 = sin(pkin(9));
t10 = t27 + (-pkin(6) * qJD(1) - t24) * t30;
t16 = t30 * qJD(2) + t32 * t24;
t39 = qJD(1) * t32;
t11 = pkin(6) * t39 + t16;
t35 = sin(qJ(4));
t42 = cos(qJ(4));
t6 = t35 * t10 + t42 * t11;
t43 = pkin(1) * t36;
t41 = cos(qJ(5));
t40 = qJD(1) * t30;
t33 = cos(pkin(8));
t38 = -pkin(1) * t33 - pkin(2);
t19 = t35 * t40 - t42 * t39;
t5 = t42 * t10 - t35 * t11;
t18 = qJD(3) + (-pkin(3) * t32 + t38) * qJD(1);
t34 = sin(qJ(5));
t23 = t38 * qJD(1) + qJD(3);
t21 = (t42 * t30 + t32 * t35) * qJD(1);
t17 = qJD(5) + t19;
t15 = -t30 * t24 + t27;
t14 = t34 * qJD(4) + t41 * t21;
t12 = -t41 * qJD(4) + t34 * t21;
t7 = t19 * pkin(4) - t21 * pkin(7) + t18;
t4 = qJD(4) * pkin(7) + t6;
t3 = -qJD(4) * pkin(4) - t5;
t2 = t34 * t7 + t41 * t4;
t1 = -t34 * t4 + t41 * t7;
t8 = [0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t33 * t43, -t31 * t43, 0, qJD(2) ^ 2 / 0.2e1 + (t31 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t36, t30 ^ 2 * t29, t30 * t36 * t32, 0, t32 ^ 2 * t29, 0, 0, -t23 * t39, t23 * t40, (-t15 * t30 + t16 * t32) * qJD(1), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * qJD(4), t19 ^ 2 / 0.2e1, -t19 * qJD(4), qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) + t18 * t19, -t6 * qJD(4) + t18 * t21, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t17, t12 ^ 2 / 0.2e1, -t12 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t12, t3 * t14 - t2 * t17, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
