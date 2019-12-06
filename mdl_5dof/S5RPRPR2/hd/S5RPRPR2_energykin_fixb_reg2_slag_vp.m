% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR2
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:36
% EndTime: 2019-12-05 17:49:36
% DurationCPUTime: 0.13s
% Computational Cost: add. (209->35), mult. (398->88), div. (0->0), fcn. (215->8), ass. (0->35)
t24 = qJD(1) + qJD(3);
t23 = t24 ^ 2;
t42 = t23 / 0.2e1;
t34 = qJD(1) ^ 2;
t41 = pkin(1) * t34;
t30 = cos(pkin(8));
t18 = (pkin(1) * t30 + pkin(2)) * qJD(1);
t32 = sin(qJ(3));
t33 = cos(qJ(3));
t28 = sin(pkin(8));
t37 = pkin(1) * qJD(1) * t28;
t12 = t32 * t18 + t33 * t37;
t10 = t24 * qJ(4) + t12;
t27 = sin(pkin(9));
t29 = cos(pkin(9));
t6 = t27 * qJD(2) + t29 * t10;
t40 = cos(qJ(5));
t39 = t24 * t27;
t38 = t24 * t29;
t11 = t33 * t18 - t32 * t37;
t36 = qJD(4) - t11;
t31 = sin(qJ(5));
t26 = t34 / 0.2e1;
t25 = qJD(2) ^ 2 / 0.2e1;
t22 = t29 * qJD(2);
t15 = (t40 * t27 + t29 * t31) * t24;
t13 = t31 * t39 - t40 * t38;
t9 = -t24 * pkin(3) + t36;
t7 = (-pkin(4) * t29 - pkin(3)) * t24 + t36;
t5 = -t27 * t10 + t22;
t4 = pkin(7) * t38 + t6;
t3 = t22 + (-pkin(7) * t24 - t10) * t27;
t2 = t31 * t3 + t40 * t4;
t1 = t40 * t3 - t31 * t4;
t8 = [0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t30 * t41, -t28 * t41, 0, t25 + (t28 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t34, 0, 0, 0, 0, 0, t42, t11 * t24, -t12 * t24, 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t25, t27 ^ 2 * t42, t27 * t23 * t29, 0, t29 ^ 2 * t42, 0, 0, -t9 * t38, t9 * t39, (-t27 * t5 + t29 * t6) * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * qJD(5), t13 ^ 2 / 0.2e1, -t13 * qJD(5), qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t7 * t13, -t2 * qJD(5) + t7 * t15, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
