% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:44
% EndTime: 2019-12-31 17:35:44
% DurationCPUTime: 0.09s
% Computational Cost: add. (74->20), mult. (151->59), div. (0->0), fcn. (78->6), ass. (0->23)
t11 = qJD(3) + qJD(4);
t10 = t11 ^ 2;
t24 = t10 / 0.2e1;
t14 = sin(qJ(4));
t17 = cos(qJ(4));
t15 = sin(qJ(3));
t22 = qJD(2) * t15;
t18 = cos(qJ(3));
t8 = qJD(3) * pkin(3) + t18 * qJD(2);
t6 = t14 * t8 + t17 * t22;
t5 = -t14 * t22 + t17 * t8;
t3 = -t11 * pkin(4) - t5;
t23 = t11 * t3;
t21 = qJD(5) * t11;
t20 = qJD(2) * qJD(3);
t19 = qJD(2) ^ 2;
t16 = cos(qJ(5));
t13 = sin(qJ(5));
t12 = qJD(1) ^ 2 / 0.2e1;
t4 = t11 * pkin(7) + t6;
t2 = -t13 * qJD(1) + t16 * t4;
t1 = -t16 * qJD(1) - t13 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 + t19 / 0.2e1, 0, 0, 0, 0, 0, qJD(3) ^ 2 / 0.2e1, t18 * t20, -t15 * t20, 0, t12 + (t15 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1) * t19, 0, 0, 0, 0, 0, t24, t5 * t11, -t6 * t11, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12, t13 ^ 2 * t24, t13 * t10 * t16, t13 * t21, t16 ^ 2 * t24, t16 * t21, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t16 * t23, -t2 * qJD(5) + t13 * t23, (-t1 * t13 + t16 * t2) * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
