% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRRR2
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:47
% EndTime: 2019-12-05 15:14:47
% DurationCPUTime: 0.15s
% Computational Cost: add. (127->34), mult. (336->90), div. (0->0), fcn. (227->8), ass. (0->31)
t31 = qJD(3) ^ 2;
t39 = t31 / 0.2e1;
t24 = sin(pkin(9));
t25 = cos(pkin(9));
t36 = qJD(1) * cos(qJ(3));
t37 = qJD(1) * sin(qJ(3));
t12 = t24 * t36 + t25 * t37;
t10 = qJD(3) * pkin(6) + t12;
t27 = sin(qJ(4));
t29 = cos(qJ(4));
t6 = t27 * qJD(2) + t29 * t10;
t38 = cos(qJ(5));
t35 = qJD(3) * t27;
t34 = qJD(3) * t29;
t33 = qJD(3) * qJD(4);
t11 = -t24 * t37 + t25 * t36;
t32 = qJD(1) ^ 2;
t26 = sin(qJ(5));
t23 = qJD(2) ^ 2 / 0.2e1;
t22 = qJD(4) + qJD(5);
t21 = t29 * qJD(2);
t15 = (t26 * t29 + t38 * t27) * qJD(3);
t13 = t26 * t35 - t38 * t34;
t9 = -qJD(3) * pkin(3) - t11;
t7 = (-pkin(4) * t29 - pkin(3)) * qJD(3) - t11;
t5 = -t27 * t10 + t21;
t4 = pkin(7) * t34 + t6;
t3 = qJD(4) * pkin(4) + t21 + (-pkin(7) * qJD(3) - t10) * t27;
t2 = t26 * t3 + t38 * t4;
t1 = -t26 * t4 + t38 * t3;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t32 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 + (t24 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1) * t32, 0, 0, 0, 0, 0, t39, t11 * qJD(3), -t12 * qJD(3), 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t23, t27 ^ 2 * t39, t27 * t31 * t29, t27 * t33, t29 ^ 2 * t39, t29 * t33, qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) - t9 * t34, -t6 * qJD(4) + t9 * t35, (-t27 * t5 + t29 * t6) * qJD(3), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t22, t13 ^ 2 / 0.2e1, -t13 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t7 * t13, t7 * t15 - t2 * t22, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
