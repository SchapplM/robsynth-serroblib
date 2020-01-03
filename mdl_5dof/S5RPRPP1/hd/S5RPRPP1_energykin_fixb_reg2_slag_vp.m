% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:13
% EndTime: 2019-12-31 18:09:13
% DurationCPUTime: 0.13s
% Computational Cost: add. (188->39), mult. (488->95), div. (0->0), fcn. (271->6), ass. (0->36)
t31 = qJD(1) ^ 2;
t25 = t31 / 0.2e1;
t26 = sin(pkin(8));
t40 = cos(pkin(8));
t27 = sin(pkin(7));
t19 = (pkin(1) * t27 + pkin(6)) * qJD(1);
t30 = cos(qJ(3));
t23 = t30 * qJD(2);
t29 = sin(qJ(3));
t35 = qJ(4) * qJD(1);
t8 = qJD(3) * pkin(3) + t23 + (-t19 - t35) * t29;
t11 = t29 * qJD(2) + t30 * t19;
t9 = t30 * t35 + t11;
t4 = t26 * t8 + t40 * t9;
t42 = pkin(1) * t31;
t38 = qJD(1) * t30;
t39 = qJD(1) * t29;
t15 = t26 * t39 - t40 * t38;
t17 = (t26 * t30 + t40 * t29) * qJD(1);
t41 = t17 * t15;
t37 = qJD(3) * t15;
t36 = t15 ^ 2 / 0.2e1;
t34 = qJD(1) * qJD(3);
t28 = cos(pkin(7));
t33 = -pkin(1) * t28 - pkin(2);
t3 = -t26 * t9 + t40 * t8;
t14 = qJD(4) + (-pkin(3) * t30 + t33) * qJD(1);
t24 = qJD(3) ^ 2 / 0.2e1;
t20 = t33 * qJD(1);
t13 = t17 * qJD(3);
t12 = t17 ^ 2 / 0.2e1;
t10 = -t29 * t19 + t23;
t5 = t15 * pkin(4) - t17 * qJ(5) + t14;
t2 = qJD(3) * qJ(5) + t4;
t1 = -qJD(3) * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t28 * t42, -t27 * t42, 0, qJD(2) ^ 2 / 0.2e1 + (t27 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t31, t29 ^ 2 * t25, t29 * t31 * t30, t29 * t34, t30 ^ 2 * t25, t30 * t34, t24, t10 * qJD(3) - t20 * t38, -t11 * qJD(3) + t20 * t39, (-t10 * t29 + t11 * t30) * qJD(1), t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t12, -t41, t13, t36, -t37, t24, t3 * qJD(3) + t14 * t15, -t4 * qJD(3) + t14 * t17, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t12, t13, t41, t24, t37, t36, -t1 * qJD(3) + t5 * t15, t1 * t17 - t2 * t15, t2 * qJD(3) - t5 * t17, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
