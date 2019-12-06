% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:15
% EndTime: 2019-12-05 15:11:15
% DurationCPUTime: 0.11s
% Computational Cost: add. (84->29), mult. (228->66), div. (0->0), fcn. (131->6), ass. (0->31)
t25 = qJD(3) ^ 2;
t35 = t25 / 0.2e1;
t34 = qJD(1) ^ 2 / 0.2e1;
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t19 = sin(pkin(8));
t33 = qJD(1) * t19;
t10 = t22 * qJD(2) + t24 * t33;
t20 = cos(pkin(8));
t32 = qJD(1) * t20;
t21 = sin(qJ(4));
t31 = qJD(3) * t21;
t23 = cos(qJ(4));
t30 = qJD(3) * t23;
t29 = qJD(3) * qJD(4);
t28 = t21 * t25 * t23;
t27 = t23 * t29;
t9 = t24 * qJD(2) - t22 * t33;
t8 = qJD(3) * pkin(6) + t10;
t4 = -t21 * t32 + t23 * t8;
t3 = -t21 * t8 - t23 * t32;
t18 = qJD(4) ^ 2 / 0.2e1;
t16 = t21 * t29;
t15 = t23 ^ 2 * t35;
t14 = t21 ^ 2 * t35;
t13 = t20 ^ 2 * t34;
t7 = -qJD(3) * pkin(3) - t9;
t5 = (-pkin(4) * t23 - qJ(5) * t21 - pkin(3)) * qJD(3) - t9;
t2 = qJD(4) * qJ(5) + t4;
t1 = -qJD(4) * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 ^ 2 * t34 + t13 + qJD(2) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t35, t9 * qJD(3), -t10 * qJD(3), 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t13, t14, t28, t16, t15, t27, t18, t3 * qJD(4) - t7 * t30, -t4 * qJD(4) + t7 * t31, (-t21 * t3 + t23 * t4) * qJD(3), t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t14, t16, -t28, t18, -t27, t15, -t1 * qJD(4) - t5 * t30, (t1 * t21 + t2 * t23) * qJD(3), t2 * qJD(4) - t5 * t31, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
