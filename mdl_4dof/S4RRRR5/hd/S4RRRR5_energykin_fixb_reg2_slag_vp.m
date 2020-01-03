% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:12
% EndTime: 2019-12-31 17:28:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (196->37), mult. (499->98), div. (0->0), fcn. (293->6), ass. (0->33)
t33 = qJD(1) ^ 2;
t43 = t33 / 0.2e1;
t42 = cos(qJ(3));
t41 = cos(qJ(4));
t32 = cos(qJ(2));
t40 = t32 * t33;
t31 = sin(qJ(2));
t14 = (-pkin(2) * t32 - pkin(6) * t31 - pkin(1)) * qJD(1);
t38 = t32 * qJD(1);
t20 = pkin(5) * t38 + qJD(2) * pkin(6);
t30 = sin(qJ(3));
t9 = t30 * t14 + t42 * t20;
t39 = qJD(1) * t31;
t37 = qJD(1) * qJD(2);
t36 = t31 * t37;
t35 = t32 * t37;
t8 = t42 * t14 - t30 * t20;
t19 = -qJD(2) * pkin(2) + pkin(5) * t39;
t23 = -qJD(3) + t38;
t29 = sin(qJ(4));
t28 = t32 ^ 2;
t27 = t31 ^ 2;
t21 = -qJD(4) + t23;
t17 = t30 * qJD(2) + t42 * t39;
t15 = -t42 * qJD(2) + t30 * t39;
t10 = t15 * pkin(3) + t19;
t7 = -t29 * t15 + t41 * t17;
t5 = t41 * t15 + t29 * t17;
t4 = -t15 * pkin(7) + t9;
t3 = -t23 * pkin(3) - t17 * pkin(7) + t8;
t2 = t29 * t3 + t41 * t4;
t1 = -t29 * t4 + t41 * t3;
t6 = [0, 0, 0, 0, 0, t43, 0, 0, 0, 0, t27 * t43, t31 * t40, t36, t28 * t43, t35, qJD(2) ^ 2 / 0.2e1, pkin(1) * t40 - pkin(5) * t36, -t33 * pkin(1) * t31 - pkin(5) * t35, (t27 + t28) * t33 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t28 / 0.2e1 + t27 / 0.2e1) * pkin(5) ^ 2) * t33, t17 ^ 2 / 0.2e1, -t17 * t15, -t17 * t23, t15 ^ 2 / 0.2e1, t15 * t23, t23 ^ 2 / 0.2e1, t19 * t15 - t8 * t23, t19 * t17 + t9 * t23, -t9 * t15 - t8 * t17, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t7 ^ 2 / 0.2e1, -t7 * t5, -t7 * t21, t5 ^ 2 / 0.2e1, t5 * t21, t21 ^ 2 / 0.2e1, -t1 * t21 + t10 * t5, t10 * t7 + t2 * t21, -t1 * t7 - t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t6;
