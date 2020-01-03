% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:46
% EndTime: 2019-12-31 19:41:47
% DurationCPUTime: 0.15s
% Computational Cost: add. (163->45), mult. (402->100), div. (0->0), fcn. (162->4), ass. (0->37)
t35 = qJD(1) ^ 2;
t45 = t35 / 0.2e1;
t44 = -pkin(2) - pkin(3);
t34 = cos(qJ(2));
t43 = t34 * t35;
t42 = qJD(1) * t34;
t14 = pkin(6) * t42 + qJD(2) * qJ(3);
t32 = sin(qJ(2));
t41 = t32 * qJD(1);
t40 = pkin(6) * t41 + qJD(3);
t39 = qJ(4) * qJD(1);
t38 = qJD(1) * qJD(2);
t9 = -qJD(1) * pkin(1) - pkin(2) * t42 - qJ(3) * t41;
t19 = t32 * t38;
t20 = t34 * t38;
t5 = pkin(3) * t42 + qJD(4) - t9;
t8 = t34 * t39 - t14;
t37 = -t32 * t39 + t40;
t33 = cos(qJ(5));
t31 = sin(qJ(5));
t29 = t34 ^ 2;
t28 = t32 ^ 2;
t26 = qJD(2) ^ 2 / 0.2e1;
t18 = t29 * t45;
t17 = t28 * t45;
t16 = qJD(5) + t41;
t15 = t32 * t43;
t13 = -qJD(2) * pkin(2) + t40;
t11 = t31 * qJD(2) + t33 * t42;
t10 = -t33 * qJD(2) + t31 * t42;
t7 = qJD(2) * pkin(4) - t8;
t6 = t44 * qJD(2) + t37;
t4 = (-pkin(7) + t44) * qJD(2) + t37;
t3 = (pkin(4) * t32 + pkin(7) * t34) * qJD(1) + t5;
t2 = t31 * t3 + t33 * t4;
t1 = t33 * t3 - t31 * t4;
t12 = [0, 0, 0, 0, 0, t45, 0, 0, 0, 0, t17, t15, t19, t18, t20, t26, pkin(1) * t43 - pkin(6) * t19, -t35 * pkin(1) * t32 - pkin(6) * t20, (t28 + t29) * t35 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t29 / 0.2e1 + t28 / 0.2e1) * pkin(6) ^ 2) * t35, t17, t19, -t15, t26, -t20, t18, -t13 * qJD(2) - t9 * t42, (t13 * t32 + t14 * t34) * qJD(1), t14 * qJD(2) - t9 * t41, t14 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t18, t15, t20, t17, t19, t26, -t8 * qJD(2) + t5 * t41, t6 * qJD(2) - t5 * t42, (-t32 * t6 + t34 * t8) * qJD(1), t6 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t10, -t11 * t16, t10 ^ 2 / 0.2e1, t10 * t16, t16 ^ 2 / 0.2e1, t1 * t16 - t7 * t10, -t7 * t11 - t2 * t16, t1 * t11 + t2 * t10, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t12;
