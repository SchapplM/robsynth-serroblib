% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:58
% EndTime: 2019-12-05 15:16:58
% DurationCPUTime: 0.12s
% Computational Cost: add. (127->34), mult. (328->88), div. (0->0), fcn. (218->8), ass. (0->33)
t29 = qJD(3) ^ 2;
t39 = t29 / 0.2e1;
t38 = qJD(1) ^ 2 / 0.2e1;
t37 = cos(qJ(5));
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t22 = sin(pkin(9));
t36 = qJD(1) * t22;
t15 = t26 * qJD(2) + t28 * t36;
t23 = cos(pkin(9));
t35 = qJD(1) * t23;
t25 = sin(qJ(4));
t34 = qJD(3) * t25;
t27 = cos(qJ(4));
t33 = qJD(3) * t27;
t32 = qJD(3) * qJD(4);
t31 = t27 * t35;
t14 = t28 * qJD(2) - t26 * t36;
t10 = qJD(3) * pkin(6) + t15;
t6 = t27 * t10 - t25 * t35;
t24 = sin(qJ(5));
t21 = qJD(4) + qJD(5);
t19 = t23 ^ 2 * t38;
t13 = (t24 * t27 + t37 * t25) * qJD(3);
t11 = t24 * t34 - t37 * t33;
t9 = -qJD(3) * pkin(3) - t14;
t7 = (-pkin(4) * t27 - pkin(3)) * qJD(3) - t14;
t5 = -t25 * t10 - t31;
t4 = pkin(7) * t33 + t6;
t3 = -t31 + qJD(4) * pkin(4) + (-pkin(7) * qJD(3) - t10) * t25;
t2 = t24 * t3 + t37 * t4;
t1 = -t24 * t4 + t37 * t3;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 ^ 2 * t38 + t19 + qJD(2) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t39, t14 * qJD(3), -t15 * qJD(3), 0, t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t19, t25 ^ 2 * t39, t25 * t29 * t27, t25 * t32, t27 ^ 2 * t39, t27 * t32, qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) - t9 * t33, -t6 * qJD(4) + t9 * t34, (-t25 * t5 + t27 * t6) * qJD(3), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t21, t11 ^ 2 / 0.2e1, -t11 * t21, t21 ^ 2 / 0.2e1, t1 * t21 + t7 * t11, t7 * t13 - t2 * t21, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
