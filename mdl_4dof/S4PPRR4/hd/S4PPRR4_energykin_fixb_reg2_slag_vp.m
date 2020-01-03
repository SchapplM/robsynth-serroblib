% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PPRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:39
% EndTime: 2019-12-31 16:18:39
% DurationCPUTime: 0.09s
% Computational Cost: add. (38->17), mult. (128->53), div. (0->0), fcn. (76->6), ass. (0->19)
t10 = sin(pkin(7));
t11 = cos(pkin(7));
t13 = sin(qJ(3));
t15 = cos(qJ(3));
t5 = (t10 * t13 - t11 * t15) * qJD(1);
t16 = qJD(3) ^ 2;
t21 = t16 / 0.2e1;
t6 = (t10 * t15 + t11 * t13) * qJD(1);
t3 = -qJD(3) * pkin(3) + t5;
t20 = qJD(3) * t3;
t19 = qJD(3) * qJD(4);
t17 = qJD(1) ^ 2;
t14 = cos(qJ(4));
t12 = sin(qJ(4));
t9 = qJD(2) ^ 2 / 0.2e1;
t4 = qJD(3) * pkin(5) + t6;
t2 = qJD(2) * t12 + t14 * t4;
t1 = qJD(2) * t14 - t12 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t17 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 + (t10 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1) * t17, 0, 0, 0, 0, 0, t21, -t5 * qJD(3), -t6 * qJD(3), 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9, t12 ^ 2 * t21, t12 * t16 * t14, t12 * t19, t14 ^ 2 * t21, t14 * t19, qJD(4) ^ 2 / 0.2e1, qJD(4) * t1 - t14 * t20, -qJD(4) * t2 + t12 * t20, (-t1 * t12 + t14 * t2) * qJD(3), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
