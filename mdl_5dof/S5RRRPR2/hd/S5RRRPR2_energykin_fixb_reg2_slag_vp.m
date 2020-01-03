% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:36
% EndTime: 2020-01-03 12:07:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (220->24), mult. (343->73), div. (0->0), fcn. (172->8), ass. (0->29)
t17 = qJD(1) + qJD(2);
t16 = qJD(3) + t17;
t15 = t16 ^ 2;
t14 = t15 / 0.2e1;
t25 = cos(qJ(2));
t31 = pkin(1) * qJD(1);
t28 = t25 * t31;
t13 = t17 * pkin(2) + t28;
t21 = sin(qJ(3));
t24 = cos(qJ(3));
t22 = sin(qJ(2));
t29 = t22 * t31;
t11 = t21 * t13 + t24 * t29;
t18 = sin(pkin(9));
t19 = cos(pkin(9));
t10 = t24 * t13 - t21 * t29;
t8 = t16 * pkin(3) + t10;
t6 = t19 * t11 + t18 * t8;
t5 = -t18 * t11 + t19 * t8;
t3 = -t16 * pkin(4) - t5;
t32 = t16 * t3;
t30 = qJD(5) * t16;
t26 = qJD(1) ^ 2;
t23 = cos(qJ(5));
t20 = sin(qJ(5));
t4 = t16 * pkin(8) + t6;
t2 = t20 * qJD(4) + t23 * t4;
t1 = t23 * qJD(4) - t20 * t4;
t7 = [0, 0, 0, 0, 0, t26 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 ^ 2 / 0.2e1, t17 * t28, -t17 * t29, 0, (t22 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t26, 0, 0, 0, 0, 0, t14, t10 * t16, -t11 * t16, 0, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t14, t5 * t16, -t6 * t16, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t20 ^ 2 * t14, t20 * t15 * t23, t20 * t30, t23 ^ 2 * t14, t23 * t30, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t23 * t32, -t2 * qJD(5) + t20 * t32, (-t1 * t20 + t2 * t23) * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
