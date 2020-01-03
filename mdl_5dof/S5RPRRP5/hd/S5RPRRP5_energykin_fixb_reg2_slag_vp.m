% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP5
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:56
% EndTime: 2019-12-31 18:40:56
% DurationCPUTime: 0.13s
% Computational Cost: add. (152->32), mult. (314->74), div. (0->0), fcn. (142->6), ass. (0->34)
t20 = qJD(1) + qJD(3);
t19 = t20 ^ 2;
t39 = t19 / 0.2e1;
t30 = qJD(1) ^ 2;
t38 = pkin(1) * t30;
t26 = sin(qJ(4));
t28 = cos(qJ(4));
t25 = cos(pkin(8));
t12 = (pkin(1) * t25 + pkin(2)) * qJD(1);
t27 = sin(qJ(3));
t29 = cos(qJ(3));
t24 = sin(pkin(8));
t33 = pkin(1) * qJD(1) * t24;
t10 = t27 * t12 + t29 * t33;
t8 = t20 * pkin(7) + t10;
t5 = t26 * qJD(2) + t28 * t8;
t37 = t20 * t26;
t36 = t20 * t28;
t35 = qJD(4) * t20;
t34 = t26 * t19 * t28;
t32 = t28 * t35;
t9 = t29 * t12 - t27 * t33;
t4 = t28 * qJD(2) - t26 * t8;
t23 = t30 / 0.2e1;
t22 = qJD(2) ^ 2 / 0.2e1;
t21 = qJD(4) ^ 2 / 0.2e1;
t17 = t26 * t35;
t16 = t28 ^ 2 * t39;
t15 = t26 ^ 2 * t39;
t7 = -t20 * pkin(3) - t9;
t3 = qJD(4) * qJ(5) + t5;
t2 = -qJD(4) * pkin(4) + qJD(5) - t4;
t1 = (-pkin(4) * t28 - qJ(5) * t26 - pkin(3)) * t20 - t9;
t6 = [0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t25 * t38, -t24 * t38, 0, t22 + (t24 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t30, 0, 0, 0, 0, 0, t39, t9 * t20, -t10 * t20, 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t22, t15, t34, t17, t16, t32, t21, t4 * qJD(4) - t7 * t36, -t5 * qJD(4) + t7 * t37, (-t26 * t4 + t28 * t5) * t20, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t15, t17, -t34, t21, -t32, t16, -t2 * qJD(4) - t1 * t36, (t2 * t26 + t28 * t3) * t20, t3 * qJD(4) - t1 * t37, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t6;
