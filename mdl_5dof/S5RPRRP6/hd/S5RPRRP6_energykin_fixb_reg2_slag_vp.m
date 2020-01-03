% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:15
% EndTime: 2019-12-31 18:43:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (192->40), mult. (466->95), div. (0->0), fcn. (249->6), ass. (0->35)
t37 = qJD(1) ^ 2;
t31 = t37 / 0.2e1;
t32 = sin(pkin(8));
t24 = (pkin(1) * t32 + pkin(6)) * qJD(1);
t35 = sin(qJ(3));
t36 = cos(qJ(3));
t16 = t35 * qJD(2) + t36 * t24;
t11 = qJD(3) * pkin(7) + t16;
t33 = cos(pkin(8));
t39 = -pkin(1) * t33 - pkin(2);
t12 = (-pkin(3) * t36 - pkin(7) * t35 + t39) * qJD(1);
t34 = sin(qJ(4));
t43 = cos(qJ(4));
t4 = t43 * t11 + t34 * t12;
t44 = pkin(1) * t37;
t42 = qJD(1) * t35;
t41 = t36 * qJD(1);
t40 = qJD(1) * qJD(3);
t3 = -t34 * t11 + t43 * t12;
t15 = t36 * qJD(2) - t35 * t24;
t10 = -qJD(3) * pkin(3) - t15;
t27 = -qJD(4) + t41;
t26 = t27 ^ 2 / 0.2e1;
t25 = t39 * qJD(1);
t23 = t34 * qJD(3) + t43 * t42;
t21 = -t43 * qJD(3) + t34 * t42;
t18 = t23 ^ 2 / 0.2e1;
t17 = t21 ^ 2 / 0.2e1;
t14 = t23 * t27;
t13 = t21 * t27;
t6 = t23 * t21;
t5 = t21 * pkin(4) + qJD(5) + t10;
t2 = -t21 * qJ(5) + t4;
t1 = -t27 * pkin(4) - t23 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t33 * t44, -t32 * t44, 0, qJD(2) ^ 2 / 0.2e1 + (t32 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t37, t35 ^ 2 * t31, t35 * t37 * t36, t35 * t40, t36 ^ 2 * t31, t36 * t40, qJD(3) ^ 2 / 0.2e1, t15 * qJD(3) - t25 * t41, -t16 * qJD(3) + t25 * t42, (-t15 * t35 + t16 * t36) * qJD(1), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t18, -t6, -t14, t17, t13, t26, t10 * t21 - t3 * t27, t10 * t23 + t4 * t27, -t4 * t21 - t3 * t23, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t18, -t6, -t14, t17, t13, t26, -t1 * t27 + t5 * t21, t2 * t27 + t5 * t23, -t1 * t23 - t2 * t21, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
