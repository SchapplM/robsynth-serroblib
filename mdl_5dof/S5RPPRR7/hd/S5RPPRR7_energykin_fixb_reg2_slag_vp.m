% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:52
% EndTime: 2019-12-31 17:59:52
% DurationCPUTime: 0.11s
% Computational Cost: add. (135->36), mult. (297->85), div. (0->0), fcn. (131->6), ass. (0->31)
t26 = qJD(1) ^ 2;
t20 = t26 / 0.2e1;
t36 = pkin(1) * t26;
t24 = sin(qJ(4));
t25 = cos(qJ(4));
t22 = cos(pkin(8));
t29 = -pkin(1) * t22 - pkin(2);
t9 = qJD(3) + (-pkin(6) + t29) * qJD(1);
t6 = t25 * qJD(2) + t24 * t9;
t35 = cos(qJ(5));
t34 = qJD(1) * t25;
t21 = sin(pkin(8));
t28 = -pkin(1) * t21 - qJ(3);
t14 = t28 * qJD(1);
t33 = t14 * qJD(1);
t32 = t24 * qJD(1);
t31 = t14 ^ 2 / 0.2e1;
t30 = qJD(1) * qJD(4);
t5 = -t24 * qJD(2) + t25 * t9;
t23 = sin(qJ(5));
t19 = qJD(2) ^ 2 / 0.2e1;
t16 = qJD(5) + t32;
t13 = t29 * qJD(1) + qJD(3);
t12 = t23 * qJD(4) + t35 * t34;
t10 = -t35 * qJD(4) + t23 * t34;
t7 = (pkin(4) * t24 - pkin(7) * t25 - t28) * qJD(1);
t4 = qJD(4) * pkin(7) + t6;
t3 = -qJD(4) * pkin(4) - t5;
t2 = t23 * t7 + t35 * t4;
t1 = -t23 * t4 + t35 * t7;
t8 = [0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t22 * t36, -t21 * t36, 0, t19 + (t21 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t26, t20, 0, 0, 0, 0, 0, 0, t13 * qJD(1), -t33, t19 + t31 + t13 ^ 2 / 0.2e1, t25 ^ 2 * t20, -t25 * t26 * t24, t25 * t30, t24 ^ 2 * t20, -t24 * t30, qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) - t14 * t32, -t6 * qJD(4) - t25 * t33, (-t24 * t6 - t25 * t5) * qJD(1), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t31, t12 ^ 2 / 0.2e1, -t12 * t10, t12 * t16, t10 ^ 2 / 0.2e1, -t10 * t16, t16 ^ 2 / 0.2e1, t1 * t16 + t3 * t10, t3 * t12 - t2 * t16, -t1 * t12 - t2 * t10, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
