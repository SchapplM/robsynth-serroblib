% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:34
% EndTime: 2020-01-03 12:00:34
% DurationCPUTime: 0.10s
% Computational Cost: add. (205->25), mult. (345->71), div. (0->0), fcn. (172->8), ass. (0->31)
t17 = qJD(1) + qJD(2);
t16 = qJD(4) + t17;
t15 = t16 ^ 2;
t34 = t15 / 0.2e1;
t26 = cos(qJ(2));
t32 = pkin(1) * qJD(1);
t29 = t26 * t32;
t13 = t17 * pkin(2) + t29;
t19 = sin(pkin(9));
t20 = cos(pkin(9));
t23 = sin(qJ(2));
t30 = t23 * t32;
t11 = t19 * t13 + t20 * t30;
t22 = sin(qJ(4));
t25 = cos(qJ(4));
t10 = t20 * t13 - t19 * t30;
t8 = t17 * pkin(3) + t10;
t6 = t25 * t11 + t22 * t8;
t5 = -t22 * t11 + t25 * t8;
t3 = -t16 * pkin(4) - t5;
t33 = t16 * t3;
t31 = qJD(5) * t16;
t27 = qJD(1) ^ 2;
t24 = cos(qJ(5));
t21 = sin(qJ(5));
t18 = qJD(3) ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t4 = t16 * pkin(8) + t6;
t2 = t21 * qJD(3) + t24 * t4;
t1 = t24 * qJD(3) - t21 * t4;
t7 = [0, 0, 0, 0, 0, t27 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t17 * t29, -t17 * t30, 0, (t23 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t27, 0, 0, 0, 0, 0, t14, t10 * t17, -t11 * t17, 0, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t18, 0, 0, 0, 0, 0, t34, t5 * t16, -t6 * t16, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18, t21 ^ 2 * t34, t21 * t15 * t24, t21 * t31, t24 ^ 2 * t34, t24 * t31, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t24 * t33, -t2 * qJD(5) + t21 * t33, (-t1 * t21 + t2 * t24) * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
