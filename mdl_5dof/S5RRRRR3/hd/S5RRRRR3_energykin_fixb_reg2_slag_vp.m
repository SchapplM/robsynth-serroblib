% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:56:28
% EndTime: 2019-12-05 18:56:28
% DurationCPUTime: 0.17s
% Computational Cost: add. (329->40), mult. (773->114), div. (0->0), fcn. (571->8), ass. (0->35)
t34 = qJD(1) ^ 2;
t43 = t34 / 0.2e1;
t42 = cos(qJ(4));
t41 = cos(qJ(5));
t32 = cos(qJ(2));
t40 = qJD(1) * t32;
t29 = sin(qJ(3));
t39 = qJD(2) * t29;
t31 = cos(qJ(3));
t38 = qJD(2) * t31;
t37 = qJD(1) * qJD(2);
t36 = t32 ^ 2 * t43;
t30 = sin(qJ(2));
t18 = t29 * t30 * qJD(1) - t31 * t40;
t20 = (t29 * t32 + t30 * t31) * qJD(1);
t10 = -pkin(1) * t40 + t18 * pkin(2) - t20 * pkin(5);
t25 = qJD(2) + qJD(3);
t21 = pkin(1) * t39 + t25 * pkin(5);
t28 = sin(qJ(4));
t7 = t42 * t10 - t28 * t21;
t22 = -pkin(1) * t38 - t25 * pkin(2);
t17 = qJD(4) + t18;
t33 = qJD(2) ^ 2;
t27 = sin(qJ(5));
t16 = qJD(5) + t17;
t15 = t42 * t20 + t28 * t25;
t13 = t28 * t20 - t42 * t25;
t11 = t13 * pkin(3) + t22;
t8 = t28 * t10 + t42 * t21;
t6 = -t27 * t13 + t41 * t15;
t4 = t41 * t13 + t27 * t15;
t3 = t17 * pkin(3) + t7;
t2 = t27 * t3 + t41 * t8;
t1 = -t27 * t8 + t41 * t3;
t5 = [0, 0, 0, 0, 0, t43, 0, 0, 0, 0, t30 ^ 2 * t43, t30 * t34 * t32, t30 * t37, t36, t32 * t37, t33 / 0.2e1, 0, 0, 0, 0, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t25, t18 ^ 2 / 0.2e1, -t18 * t25, t25 ^ 2 / 0.2e1, (-t18 * t40 + t25 * t38) * pkin(1), (-t20 * t40 - t25 * t39) * pkin(1), (-t18 * t29 - t20 * t31) * qJD(2) * pkin(1), (t36 + (t29 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1) * t33) * pkin(1) ^ 2, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t17, t13 ^ 2 / 0.2e1, -t13 * t17, t17 ^ 2 / 0.2e1, t22 * t13 + t7 * t17, t22 * t15 - t8 * t17, -t8 * t13 - t7 * t15, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t6 ^ 2 / 0.2e1, -t6 * t4, t6 * t16, t4 ^ 2 / 0.2e1, -t4 * t16, t16 ^ 2 / 0.2e1, t1 * t16 + t11 * t4, t11 * t6 - t2 * t16, -t1 * t6 - t2 * t4, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1;];
T_reg = t5;
