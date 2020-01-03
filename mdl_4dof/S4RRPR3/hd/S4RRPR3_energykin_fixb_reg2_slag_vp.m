% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:33
% EndTime: 2019-12-31 17:01:33
% DurationCPUTime: 0.09s
% Computational Cost: add. (88->18), mult. (177->60), div. (0->0), fcn. (78->6), ass. (0->23)
t12 = qJD(1) + qJD(2);
t11 = t12 ^ 2;
t10 = t11 / 0.2e1;
t13 = sin(pkin(7));
t14 = cos(pkin(7));
t16 = sin(qJ(2));
t24 = pkin(1) * qJD(1);
t22 = t16 * t24;
t18 = cos(qJ(2));
t21 = t18 * t24;
t8 = t12 * pkin(2) + t21;
t6 = t13 * t8 + t14 * t22;
t5 = -t13 * t22 + t14 * t8;
t3 = -t12 * pkin(3) - t5;
t25 = t12 * t3;
t23 = qJD(4) * t12;
t19 = qJD(1) ^ 2;
t17 = cos(qJ(4));
t15 = sin(qJ(4));
t4 = pkin(6) * t12 + t6;
t2 = t15 * qJD(3) + t17 * t4;
t1 = t17 * qJD(3) - t15 * t4;
t7 = [0, 0, 0, 0, 0, t19 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12 * t21, -t12 * t22, 0, (t16 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t19, 0, 0, 0, 0, 0, t10, t5 * t12, -t6 * t12, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t15 ^ 2 * t10, t15 * t11 * t17, t15 * t23, t17 ^ 2 * t10, t17 * t23, qJD(4) ^ 2 / 0.2e1, t1 * qJD(4) - t17 * t25, -qJD(4) * t2 + t15 * t25, (-t1 * t15 + t17 * t2) * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
