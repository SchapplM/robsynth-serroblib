% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:21
% EndTime: 2019-12-05 15:03:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (65->23), mult. (164->56), div. (0->0), fcn. (92->6), ass. (0->22)
t12 = sin(pkin(8));
t13 = cos(pkin(8));
t15 = sin(qJ(3));
t17 = cos(qJ(3));
t8 = (t12 * t17 + t13 * t15) * qJD(1);
t18 = qJD(3) ^ 2;
t10 = t18 / 0.2e1;
t5 = qJD(3) * qJ(4) + t8;
t24 = t5 * qJD(3);
t23 = t5 ^ 2 / 0.2e1;
t22 = qJD(3) * qJD(5);
t7 = (-t12 * t15 + t13 * t17) * qJD(1);
t20 = qJD(4) - t7;
t19 = qJD(1) ^ 2;
t16 = cos(qJ(5));
t14 = sin(qJ(5));
t11 = qJD(2) ^ 2 / 0.2e1;
t4 = -qJD(3) * pkin(3) + t20;
t3 = (-pkin(3) - pkin(6)) * qJD(3) + t20;
t2 = t16 * qJD(2) + t14 * t3;
t1 = -t14 * qJD(2) + t16 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t19 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + (t12 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1) * t19, 0, 0, 0, 0, 0, t10, t7 * qJD(3), -t8 * qJD(3), 0, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t11, t10, 0, 0, 0, 0, 0, 0, t4 * qJD(3), t24, t11 + t23 + t4 ^ 2 / 0.2e1, t16 ^ 2 * t10, -t16 * t18 * t14, t16 * t22, t14 ^ 2 * t10, -t14 * t22, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) + t14 * t24, -t2 * qJD(5) + t16 * t24, (-t1 * t16 - t14 * t2) * qJD(3), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t23;];
T_reg = t6;
