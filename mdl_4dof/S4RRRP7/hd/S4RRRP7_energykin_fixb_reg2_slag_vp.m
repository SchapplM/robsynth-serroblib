% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:03
% EndTime: 2019-12-31 17:21:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (126->33), mult. (341->78), div. (0->0), fcn. (174->4), ass. (0->32)
t26 = qJD(1) ^ 2;
t38 = t26 / 0.2e1;
t25 = cos(qJ(2));
t32 = t25 * qJD(1);
t15 = pkin(5) * t32 + qJD(2) * pkin(6);
t23 = sin(qJ(3));
t37 = cos(qJ(3));
t24 = sin(qJ(2));
t8 = (-pkin(2) * t25 - pkin(6) * t24 - pkin(1)) * qJD(1);
t5 = t37 * t15 + t23 * t8;
t33 = qJD(1) * t24;
t10 = -t37 * qJD(2) + t23 * t33;
t12 = t23 * qJD(2) + t37 * t33;
t36 = t12 * t10;
t17 = -qJD(3) + t32;
t35 = t17 * t10;
t34 = t25 * t26;
t31 = t10 ^ 2 / 0.2e1;
t30 = qJD(1) * qJD(2);
t29 = t24 * t30;
t28 = t25 * t30;
t14 = -qJD(2) * pkin(2) + pkin(5) * t33;
t4 = -t23 * t15 + t37 * t8;
t22 = t25 ^ 2;
t21 = t24 ^ 2;
t16 = t17 ^ 2 / 0.2e1;
t9 = t12 ^ 2 / 0.2e1;
t6 = t12 * t17;
t3 = t10 * pkin(3) - t12 * qJ(4) + t14;
t2 = -t17 * qJ(4) + t5;
t1 = t17 * pkin(3) + qJD(4) - t4;
t7 = [0, 0, 0, 0, 0, t38, 0, 0, 0, 0, t21 * t38, t24 * t34, t29, t22 * t38, t28, qJD(2) ^ 2 / 0.2e1, pkin(1) * t34 - pkin(5) * t29, -t26 * pkin(1) * t24 - pkin(5) * t28, (t21 + t22) * t26 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t22 / 0.2e1 + t21 / 0.2e1) * pkin(5) ^ 2) * t26, t9, -t36, -t6, t31, t35, t16, t14 * t10 - t4 * t17, t14 * t12 + t5 * t17, -t5 * t10 - t4 * t12, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t9, -t6, t36, t16, -t35, t31, t1 * t17 + t3 * t10, t1 * t12 - t2 * t10, -t3 * t12 - t2 * t17, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t7;
