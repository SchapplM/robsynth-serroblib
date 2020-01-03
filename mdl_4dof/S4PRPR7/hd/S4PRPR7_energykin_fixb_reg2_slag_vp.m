% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:54
% EndTime: 2019-12-31 16:25:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (38->15), mult. (97->44), div. (0->0), fcn. (31->4), ass. (0->19)
t13 = qJD(2) ^ 2;
t6 = t13 / 0.2e1;
t12 = cos(qJ(2));
t15 = -t12 * qJD(1) + qJD(3);
t2 = (-pkin(2) - pkin(5)) * qJD(2) + t15;
t20 = qJD(4) * t2;
t10 = sin(qJ(2));
t4 = qJD(2) * qJ(3) + qJD(1) * t10;
t19 = t4 * qJD(2);
t18 = t4 ^ 2 / 0.2e1;
t17 = qJD(1) * qJD(2);
t16 = qJD(2) * qJD(4);
t14 = qJD(1) ^ 2;
t11 = cos(qJ(4));
t9 = sin(qJ(4));
t8 = t11 ^ 2;
t7 = t9 ^ 2;
t3 = -qJD(2) * pkin(2) + t15;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t14 / 0.2e1, 0, 0, 0, 0, 0, t6, t12 * t17, -t10 * t17, 0, (t10 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1) * t14, t6, 0, 0, 0, 0, 0, 0, t3 * qJD(2), t19, t18 + t3 ^ 2 / 0.2e1, t8 * t6, -t11 * t13 * t9, t11 * t16, t7 * t6, -t9 * t16, qJD(4) ^ 2 / 0.2e1, t11 * t20 + t19 * t9, t11 * t19 - t20 * t9, (-t7 - t8) * t2 * qJD(2), t18 + (t7 / 0.2e1 + t8 / 0.2e1) * t2 ^ 2;];
T_reg = t1;
