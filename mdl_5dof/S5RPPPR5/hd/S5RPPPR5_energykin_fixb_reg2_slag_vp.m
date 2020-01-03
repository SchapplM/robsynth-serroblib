% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:29
% EndTime: 2019-12-31 17:46:29
% DurationCPUTime: 0.11s
% Computational Cost: add. (169->32), mult. (337->83), div. (0->0), fcn. (159->6), ass. (0->27)
t31 = qJD(1) ^ 2;
t23 = t31 / 0.2e1;
t17 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t25 = sin(pkin(7));
t27 = cos(pkin(7));
t32 = qJ(2) * qJD(1);
t12 = t25 * t17 + t27 * t32;
t10 = -qJD(1) * qJ(4) + t12;
t24 = sin(pkin(8));
t26 = cos(pkin(8));
t6 = t24 * qJD(3) + t26 * t10;
t33 = qJD(1) * t26;
t11 = t27 * t17 - t25 * t32;
t9 = qJD(1) * pkin(3) + qJD(4) - t11;
t30 = cos(qJ(5));
t29 = sin(qJ(5));
t22 = t26 * qJD(3);
t20 = -pkin(1) * qJD(1) + qJD(2);
t14 = (-t24 * t30 - t26 * t29) * qJD(1);
t13 = (t24 * t29 - t26 * t30) * qJD(1);
t7 = pkin(4) * t33 + t9;
t5 = -t24 * t10 + t22;
t4 = -pkin(6) * t33 + t6;
t3 = t22 + (pkin(6) * qJD(1) - t10) * t24;
t2 = t29 * t3 + t30 * t4;
t1 = -t29 * t4 + t30 * t3;
t8 = [0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, -t20 * qJD(1), 0, t31 * qJ(2), qJ(2) ^ 2 * t23 + t20 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t23, -t11 * qJD(1), t12 * qJD(1), 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t24 ^ 2 * t23, t24 * t31 * t26, 0, t26 ^ 2 * t23, 0, 0, t9 * t33, -t9 * t24 * qJD(1), (t24 * t5 - t26 * t6) * qJD(1), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, t14 * t13, t14 * qJD(5), t13 ^ 2 / 0.2e1, t13 * qJD(5), qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t7 * t13, -t2 * qJD(5) + t7 * t14, -t1 * t14 + t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
