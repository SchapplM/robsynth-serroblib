% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:13
% EndTime: 2019-12-05 15:07:13
% DurationCPUTime: 0.13s
% Computational Cost: add. (82->28), mult. (239->69), div. (0->0), fcn. (143->6), ass. (0->31)
t29 = qJD(3) ^ 2;
t37 = t29 / 0.2e1;
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t35 = qJD(1) * cos(qJ(3));
t36 = qJD(1) * sin(qJ(3));
t10 = t23 * t35 + t24 * t36;
t8 = qJD(3) * pkin(6) + t10;
t4 = t25 * qJD(2) + t27 * t8;
t34 = qJD(3) * t25;
t33 = qJD(3) * t27;
t32 = qJ(5) * qJD(3);
t31 = qJD(3) * qJD(4);
t9 = -t23 * t36 + t24 * t35;
t30 = qJD(1) ^ 2;
t22 = qJD(2) ^ 2 / 0.2e1;
t21 = qJD(4) ^ 2 / 0.2e1;
t20 = t27 * qJD(2);
t18 = t27 * t31;
t17 = t25 * t31;
t16 = t27 ^ 2 * t37;
t15 = t25 ^ 2 * t37;
t14 = t25 * t29 * t27;
t7 = -qJD(3) * pkin(3) - t9;
t5 = qJD(5) + (-pkin(4) * t27 - pkin(3)) * qJD(3) - t9;
t3 = -t25 * t8 + t20;
t2 = t27 * t32 + t4;
t1 = qJD(4) * pkin(4) + t20 + (-t8 - t32) * t25;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t30 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 + (t23 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1) * t30, 0, 0, 0, 0, 0, t37, t9 * qJD(3), -t10 * qJD(3), 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t22, t15, t14, t17, t16, t18, t21, t3 * qJD(4) - t7 * t33, -t4 * qJD(4) + t7 * t34, (-t25 * t3 + t27 * t4) * qJD(3), t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t15, t14, t17, t16, t18, t21, t1 * qJD(4) - t5 * t33, -t2 * qJD(4) + t5 * t34, (-t1 * t25 + t2 * t27) * qJD(3), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
