% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:48
% EndTime: 2019-12-31 16:39:48
% DurationCPUTime: 0.08s
% Computational Cost: add. (68->22), mult. (144->55), div. (0->0), fcn. (48->4), ass. (0->19)
t16 = qJD(1) ^ 2;
t11 = t16 / 0.2e1;
t12 = sin(pkin(6));
t13 = cos(pkin(6));
t18 = qJ(2) * qJD(1);
t8 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t6 = t12 * t8 + t13 * t18;
t20 = t13 * t8;
t3 = -t20 + (qJ(2) * t12 + pkin(3)) * qJD(1);
t19 = qJD(1) * t3;
t17 = qJD(1) * qJD(4);
t15 = cos(qJ(4));
t14 = sin(qJ(4));
t10 = -qJD(1) * pkin(1) + qJD(2);
t5 = -t12 * t18 + t20;
t4 = -qJD(1) * pkin(5) + t6;
t2 = qJD(3) * t14 + t15 * t4;
t1 = qJD(3) * t15 - t14 * t4;
t7 = [0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, -t10 * qJD(1), 0, t16 * qJ(2), qJ(2) ^ 2 * t11 + t10 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t11, -t5 * qJD(1), t6 * qJD(1), 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t14 ^ 2 * t11, t14 * t16 * t15, -t14 * t17, t15 ^ 2 * t11, -t15 * t17, qJD(4) ^ 2 / 0.2e1, qJD(4) * t1 + t15 * t19, -qJD(4) * t2 - t14 * t19, (t1 * t14 - t15 * t2) * qJD(1), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
