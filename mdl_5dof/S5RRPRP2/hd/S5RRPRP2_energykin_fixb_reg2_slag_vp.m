% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:42
% EndTime: 2019-12-31 19:49:42
% DurationCPUTime: 0.14s
% Computational Cost: add. (176->31), mult. (312->75), div. (0->0), fcn. (142->6), ass. (0->33)
t21 = qJD(1) + qJD(2);
t20 = t21 ^ 2;
t18 = t20 / 0.2e1;
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t28 = cos(qJ(2));
t36 = pkin(1) * qJD(1);
t32 = t28 * t36;
t12 = t21 * pkin(2) + t32;
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t26 = sin(qJ(2));
t33 = t26 * t36;
t10 = t23 * t12 + t24 * t33;
t8 = t21 * pkin(7) + t10;
t5 = t25 * qJD(3) + t27 * t8;
t38 = t21 * t25;
t37 = t21 * t27;
t35 = qJD(4) * t21;
t34 = t25 * t20 * t27;
t31 = t27 * t35;
t9 = t24 * t12 - t23 * t33;
t4 = t27 * qJD(3) - t25 * t8;
t29 = qJD(1) ^ 2;
t22 = qJD(4) ^ 2 / 0.2e1;
t17 = t25 * t35;
t16 = t27 ^ 2 * t18;
t15 = t25 ^ 2 * t18;
t7 = -t21 * pkin(3) - t9;
t3 = qJD(4) * qJ(5) + t5;
t2 = -qJD(4) * pkin(4) + qJD(5) - t4;
t1 = (-pkin(4) * t27 - qJ(5) * t25 - pkin(3)) * t21 - t9;
t6 = [0, 0, 0, 0, 0, t29 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t21 * t32, -t21 * t33, 0, (t26 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t29, 0, 0, 0, 0, 0, t18, t9 * t21, -t10 * t21, 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t15, t34, t17, t16, t31, t22, t4 * qJD(4) - t7 * t37, -t5 * qJD(4) + t7 * t38, (-t25 * t4 + t27 * t5) * t21, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t15, t17, -t34, t22, -t31, t16, -t2 * qJD(4) - t1 * t37, (t2 * t25 + t27 * t3) * t21, t3 * qJD(4) - t1 * t38, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t6;
