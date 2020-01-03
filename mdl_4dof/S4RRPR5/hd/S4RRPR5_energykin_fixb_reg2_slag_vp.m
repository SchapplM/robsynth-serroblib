% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:32
% EndTime: 2019-12-31 17:03:32
% DurationCPUTime: 0.08s
% Computational Cost: add. (72->16), mult. (116->46), div. (0->0), fcn. (31->4), ass. (0->24)
t11 = sin(qJ(4));
t9 = t11 ^ 2;
t26 = t9 / 0.2e1;
t12 = sin(qJ(2));
t23 = pkin(1) * qJD(1);
t19 = t12 * t23;
t8 = qJD(1) + qJD(2);
t4 = t8 * qJ(3) + t19;
t25 = t4 * t8;
t13 = cos(qJ(4));
t10 = t13 ^ 2;
t24 = t10 / 0.2e1;
t22 = t4 ^ 2 / 0.2e1;
t21 = qJD(4) * t11;
t20 = qJD(4) * t13;
t14 = cos(qJ(2));
t18 = t14 * t23;
t17 = qJD(3) - t18;
t15 = qJD(1) ^ 2;
t7 = t8 ^ 2;
t6 = t7 / 0.2e1;
t3 = -pkin(2) * t8 + t17;
t2 = (-pkin(2) - pkin(6)) * t8 + t17;
t1 = [0, 0, 0, 0, 0, t15 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t8 * t18, -t8 * t19, 0, (t12 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t15, t6, 0, 0, 0, 0, 0, 0, t3 * t8, t25, t22 + t3 ^ 2 / 0.2e1, t7 * t24, -t13 * t7 * t11, t8 * t20, t7 * t26, -t8 * t21, qJD(4) ^ 2 / 0.2e1, t11 * t25 + t2 * t20, t13 * t25 - t2 * t21, (-t10 - t9) * t8 * t2, t22 + (t26 + t24) * t2 ^ 2;];
T_reg = t1;
