% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRPR4
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:59
% EndTime: 2019-12-31 16:21:59
% DurationCPUTime: 0.07s
% Computational Cost: add. (29->14), mult. (82->36), div. (0->0), fcn. (22->2), ass. (0->13)
t10 = qJD(2) ^ 2;
t6 = t10 / 0.2e1;
t12 = t10 * qJ(3);
t11 = qJD(2) * qJD(4);
t9 = cos(qJ(4));
t8 = sin(qJ(4));
t7 = qJD(1) ^ 2 / 0.2e1;
t5 = qJ(3) ^ 2 * t6;
t4 = -qJD(2) * pkin(2) + qJD(3);
t3 = qJD(3) + (-pkin(2) - pkin(5)) * qJD(2);
t2 = t9 * qJD(1) + t8 * t3;
t1 = -t8 * qJD(1) + t9 * t3;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, t6, 0, 0, 0, t7, t6, 0, 0, 0, 0, 0, 0, t4 * qJD(2), t12, t7 + t5 + t4 ^ 2 / 0.2e1, t9 ^ 2 * t6, -t9 * t10 * t8, t9 * t11, t8 ^ 2 * t6, -t8 * t11, qJD(4) ^ 2 / 0.2e1, t1 * qJD(4) + t8 * t12, -t2 * qJD(4) + t9 * t12, (-t1 * t9 - t2 * t8) * qJD(2), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5;];
T_reg = t13;
