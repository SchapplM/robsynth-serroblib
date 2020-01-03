% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:45
% EndTime: 2019-12-31 16:33:45
% DurationCPUTime: 0.08s
% Computational Cost: add. (65->15), mult. (139->54), div. (0->0), fcn. (70->6), ass. (0->24)
t13 = sin(qJ(4));
t11 = t13 ^ 2;
t26 = t11 / 0.2e1;
t16 = cos(qJ(4));
t12 = t16 ^ 2;
t25 = t12 / 0.2e1;
t14 = sin(qJ(3));
t17 = cos(qJ(3));
t15 = sin(qJ(2));
t23 = qJD(1) * t15;
t18 = cos(qJ(2));
t7 = qJD(2) * pkin(2) + qJD(1) * t18;
t5 = t14 * t7 + t17 * t23;
t10 = qJD(2) + qJD(3);
t4 = -t14 * t23 + t17 * t7;
t2 = -pkin(3) * t10 - t4;
t24 = t10 * t2;
t22 = qJD(4) * t13;
t21 = qJD(4) * t16;
t20 = qJD(1) * qJD(2);
t19 = qJD(1) ^ 2;
t9 = t10 ^ 2;
t3 = pkin(6) * t10 + t5;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t19 / 0.2e1, 0, 0, 0, 0, 0, qJD(2) ^ 2 / 0.2e1, t18 * t20, -t15 * t20, 0, (t15 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1) * t19, 0, 0, 0, 0, 0, t9 / 0.2e1, t4 * t10, -t5 * t10, 0, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t9 * t26, t13 * t9 * t16, t10 * t22, t9 * t25, t10 * t21, qJD(4) ^ 2 / 0.2e1, -t16 * t24 - t3 * t22, t13 * t24 - t3 * t21, (t11 + t12) * t3 * t10, t2 ^ 2 / 0.2e1 + (t25 + t26) * t3 ^ 2;];
T_reg = t1;
