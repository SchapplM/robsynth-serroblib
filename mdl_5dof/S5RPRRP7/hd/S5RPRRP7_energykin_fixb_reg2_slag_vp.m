% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:35
% EndTime: 2019-12-31 18:45:36
% DurationCPUTime: 0.13s
% Computational Cost: add. (192->40), mult. (460->95), div. (0->0), fcn. (243->6), ass. (0->35)
t33 = qJD(1) ^ 2;
t27 = t33 / 0.2e1;
t31 = sin(qJ(3));
t32 = cos(qJ(3));
t29 = cos(pkin(8));
t35 = -pkin(1) * t29 - pkin(2);
t10 = (-pkin(3) * t32 - pkin(7) * t31 + t35) * qJD(1);
t30 = sin(qJ(4));
t42 = cos(qJ(4));
t28 = sin(pkin(8));
t20 = (pkin(1) * t28 + pkin(6)) * qJD(1);
t13 = t31 * qJD(2) + t32 * t20;
t9 = qJD(3) * pkin(7) + t13;
t5 = t30 * t10 + t42 * t9;
t43 = pkin(1) * t33;
t39 = qJD(1) * t31;
t17 = -t42 * qJD(3) + t30 * t39;
t19 = t30 * qJD(3) + t42 * t39;
t41 = t19 * t17;
t38 = t32 * qJD(1);
t23 = -qJD(4) + t38;
t40 = t23 * t17;
t37 = t17 ^ 2 / 0.2e1;
t36 = qJD(1) * qJD(3);
t12 = t32 * qJD(2) - t31 * t20;
t4 = t42 * t10 - t30 * t9;
t8 = -qJD(3) * pkin(3) - t12;
t22 = t23 ^ 2 / 0.2e1;
t21 = t35 * qJD(1);
t14 = t19 ^ 2 / 0.2e1;
t11 = t19 * t23;
t3 = t17 * pkin(4) - t19 * qJ(5) + t8;
t2 = -t23 * qJ(5) + t5;
t1 = t23 * pkin(4) + qJD(5) - t4;
t6 = [0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29 * t43, -t28 * t43, 0, qJD(2) ^ 2 / 0.2e1 + (t28 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t33, t31 ^ 2 * t27, t31 * t33 * t32, t31 * t36, t32 ^ 2 * t27, t32 * t36, qJD(3) ^ 2 / 0.2e1, t12 * qJD(3) - t21 * t38, -t13 * qJD(3) + t21 * t39, (-t12 * t31 + t13 * t32) * qJD(1), t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14, -t41, -t11, t37, t40, t22, t8 * t17 - t4 * t23, t8 * t19 + t5 * t23, -t5 * t17 - t4 * t19, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t14, -t11, t41, t22, -t40, t37, t1 * t23 + t3 * t17, t1 * t19 - t2 * t17, -t3 * t19 - t2 * t23, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
