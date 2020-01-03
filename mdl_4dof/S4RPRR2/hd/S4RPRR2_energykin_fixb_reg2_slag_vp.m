% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:12
% EndTime: 2019-12-31 16:48:12
% DurationCPUTime: 0.09s
% Computational Cost: add. (73->19), mult. (179->59), div. (0->0), fcn. (78->6), ass. (0->24)
t11 = qJD(1) + qJD(3);
t10 = t11 ^ 2;
t26 = t10 / 0.2e1;
t17 = sin(qJ(3));
t19 = cos(qJ(3));
t14 = sin(pkin(7));
t22 = pkin(1) * qJD(1) * t14;
t15 = cos(pkin(7));
t8 = (pkin(1) * t15 + pkin(2)) * qJD(1);
t6 = t17 * t8 + t19 * t22;
t20 = qJD(1) ^ 2;
t25 = pkin(1) * t20;
t5 = -t17 * t22 + t19 * t8;
t3 = -t11 * pkin(3) - t5;
t24 = t11 * t3;
t23 = qJD(4) * t11;
t18 = cos(qJ(4));
t16 = sin(qJ(4));
t13 = t20 / 0.2e1;
t12 = qJD(2) ^ 2 / 0.2e1;
t4 = pkin(6) * t11 + t6;
t2 = t16 * qJD(2) + t18 * t4;
t1 = qJD(2) * t18 - t16 * t4;
t7 = [0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t15 * t25, -t14 * t25, 0, t12 + (t14 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t20, 0, 0, 0, 0, 0, t26, t5 * t11, -t6 * t11, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12, t16 ^ 2 * t26, t16 * t10 * t18, t16 * t23, t18 ^ 2 * t26, t18 * t23, qJD(4) ^ 2 / 0.2e1, t1 * qJD(4) - t18 * t24, -qJD(4) * t2 + t16 * t24, (-t1 * t16 + t18 * t2) * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
