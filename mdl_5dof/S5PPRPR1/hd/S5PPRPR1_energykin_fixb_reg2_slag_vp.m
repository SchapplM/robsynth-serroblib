% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:18
% EndTime: 2019-12-05 15:01:18
% DurationCPUTime: 0.12s
% Computational Cost: add. (115->32), mult. (308->82), div. (0->0), fcn. (213->8), ass. (0->30)
t30 = qJD(3) ^ 2;
t38 = t30 / 0.2e1;
t24 = sin(pkin(8));
t26 = cos(pkin(8));
t35 = qJD(1) * cos(qJ(3));
t36 = qJD(1) * sin(qJ(3));
t15 = t24 * t35 + t26 * t36;
t10 = qJD(3) * qJ(4) + t15;
t23 = sin(pkin(9));
t25 = cos(pkin(9));
t6 = t23 * qJD(2) + t25 * t10;
t37 = cos(qJ(5));
t34 = qJD(3) * t23;
t33 = qJD(3) * t25;
t13 = -t24 * t36 + t26 * t35;
t32 = qJD(4) - t13;
t31 = qJD(1) ^ 2;
t27 = sin(qJ(5));
t22 = qJD(2) ^ 2 / 0.2e1;
t21 = t25 * qJD(2);
t14 = (t37 * t23 + t25 * t27) * qJD(3);
t11 = t27 * t34 - t37 * t33;
t9 = -qJD(3) * pkin(3) + t32;
t7 = (-pkin(4) * t25 - pkin(3)) * qJD(3) + t32;
t5 = -t23 * t10 + t21;
t4 = pkin(6) * t33 + t6;
t3 = t21 + (-pkin(6) * qJD(3) - t10) * t23;
t2 = t27 * t3 + t37 * t4;
t1 = -t27 * t4 + t37 * t3;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t31 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 + (t24 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1) * t31, 0, 0, 0, 0, 0, t38, t13 * qJD(3), -t15 * qJD(3), 0, t15 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t22, t23 ^ 2 * t38, t23 * t30 * t25, 0, t25 ^ 2 * t38, 0, 0, -t9 * t33, t9 * t34, (-t23 * t5 + t25 * t6) * qJD(3), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t11, t14 * qJD(5), t11 ^ 2 / 0.2e1, -t11 * qJD(5), qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t7 * t11, -t2 * qJD(5) + t7 * t14, -t1 * t14 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
