% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRRP2
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:09
% EndTime: 2019-12-05 15:09:09
% DurationCPUTime: 0.11s
% Computational Cost: add. (84->29), mult. (236->68), div. (0->0), fcn. (140->6), ass. (0->29)
t26 = qJD(3) ^ 2;
t35 = t26 / 0.2e1;
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t33 = qJD(1) * cos(qJ(3));
t34 = qJD(1) * sin(qJ(3));
t10 = t20 * t33 + t21 * t34;
t8 = qJD(3) * pkin(6) + t10;
t5 = t22 * qJD(2) + t24 * t8;
t32 = qJD(3) * t22;
t31 = qJD(3) * t24;
t30 = qJD(3) * qJD(4);
t29 = t22 * t26 * t24;
t28 = t24 * t30;
t9 = -t20 * t34 + t21 * t33;
t4 = t24 * qJD(2) - t22 * t8;
t27 = qJD(1) ^ 2;
t19 = qJD(2) ^ 2 / 0.2e1;
t18 = qJD(4) ^ 2 / 0.2e1;
t16 = t22 * t30;
t15 = t24 ^ 2 * t35;
t14 = t22 ^ 2 * t35;
t7 = -qJD(3) * pkin(3) - t9;
t3 = (-pkin(4) * t24 - qJ(5) * t22 - pkin(3)) * qJD(3) - t9;
t2 = qJD(4) * qJ(5) + t5;
t1 = -qJD(4) * pkin(4) + qJD(5) - t4;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t27 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 + (t20 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1) * t27, 0, 0, 0, 0, 0, t35, t9 * qJD(3), -t10 * qJD(3), 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t19, t14, t29, t16, t15, t28, t18, t4 * qJD(4) - t7 * t31, -t5 * qJD(4) + t7 * t32, (-t22 * t4 + t24 * t5) * qJD(3), t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t14, t16, -t29, t18, -t28, t15, -t1 * qJD(4) - t3 * t31, (t1 * t22 + t2 * t24) * qJD(3), t2 * qJD(4) - t3 * t32, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
