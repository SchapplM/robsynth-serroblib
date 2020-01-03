% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:50
% EndTime: 2019-12-31 18:19:50
% DurationCPUTime: 0.14s
% Computational Cost: add. (259->45), mult. (639->114), div. (0->0), fcn. (383->8), ass. (0->38)
t37 = qJD(1) ^ 2;
t30 = t37 / 0.2e1;
t32 = sin(pkin(8));
t23 = (pkin(1) * t32 + pkin(6)) * qJD(1);
t36 = cos(qJ(3));
t28 = t36 * qJD(2);
t35 = sin(qJ(3));
t41 = qJ(4) * qJD(1);
t10 = qJD(3) * pkin(3) + t28 + (-t23 - t41) * t35;
t16 = t35 * qJD(2) + t36 * t23;
t14 = t36 * t41 + t16;
t31 = sin(pkin(9));
t44 = cos(pkin(9));
t6 = t31 * t10 + t44 * t14;
t46 = pkin(1) * t37;
t45 = cos(qJ(5));
t43 = qJD(1) * t35;
t42 = qJD(1) * t36;
t40 = qJD(1) * qJD(3);
t33 = cos(pkin(8));
t39 = -pkin(1) * t33 - pkin(2);
t19 = t31 * t43 - t44 * t42;
t5 = t44 * t10 - t31 * t14;
t18 = qJD(4) + (-pkin(3) * t36 + t39) * qJD(1);
t34 = sin(qJ(5));
t29 = qJD(3) ^ 2 / 0.2e1;
t24 = t39 * qJD(1);
t21 = (t31 * t36 + t44 * t35) * qJD(1);
t17 = qJD(5) + t19;
t15 = -t35 * t23 + t28;
t13 = t34 * qJD(3) + t45 * t21;
t11 = -t45 * qJD(3) + t34 * t21;
t7 = t19 * pkin(4) - t21 * pkin(7) + t18;
t4 = qJD(3) * pkin(7) + t6;
t3 = -qJD(3) * pkin(4) - t5;
t2 = t34 * t7 + t45 * t4;
t1 = -t34 * t4 + t45 * t7;
t8 = [0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t33 * t46, -t32 * t46, 0, qJD(2) ^ 2 / 0.2e1 + (t32 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t37, t35 ^ 2 * t30, t35 * t37 * t36, t35 * t40, t36 ^ 2 * t30, t36 * t40, t29, t15 * qJD(3) - t24 * t42, -t16 * qJD(3) + t24 * t43, (-t15 * t35 + t16 * t36) * qJD(1), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * qJD(3), t19 ^ 2 / 0.2e1, -t19 * qJD(3), t29, t5 * qJD(3) + t18 * t19, -t6 * qJD(3) + t18 * t21, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t17, t11 ^ 2 / 0.2e1, -t11 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t11, t3 * t13 - t2 * t17, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
