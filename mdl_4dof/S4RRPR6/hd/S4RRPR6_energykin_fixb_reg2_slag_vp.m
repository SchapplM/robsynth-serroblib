% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPR6
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:53
% EndTime: 2019-12-31 17:04:53
% DurationCPUTime: 0.12s
% Computational Cost: add. (186->36), mult. (546->94), div. (0->0), fcn. (340->6), ass. (0->34)
t30 = qJD(1) ^ 2;
t41 = t30 / 0.2e1;
t40 = cos(qJ(4));
t29 = cos(qJ(2));
t39 = t29 * t30;
t38 = pkin(5) + qJ(3);
t28 = sin(qJ(2));
t36 = qJD(1) * t28;
t18 = qJD(2) * pkin(2) - t38 * t36;
t35 = qJD(1) * t29;
t19 = t38 * t35;
t26 = sin(pkin(7));
t37 = cos(pkin(7));
t9 = t26 * t18 + t37 * t19;
t34 = qJD(1) * qJD(2);
t33 = t28 * t34;
t32 = t29 * t34;
t8 = t37 * t18 - t26 * t19;
t20 = qJD(3) + (-pkin(2) * t29 - pkin(1)) * qJD(1);
t27 = sin(qJ(4));
t25 = t29 ^ 2;
t24 = t28 ^ 2;
t23 = qJD(2) ^ 2 / 0.2e1;
t22 = qJD(2) + qJD(4);
t16 = (t26 * t29 + t37 * t28) * qJD(1);
t14 = t26 * t36 - t37 * t35;
t10 = t14 * pkin(3) + t20;
t7 = -t27 * t14 + t40 * t16;
t5 = t40 * t14 + t27 * t16;
t4 = -t14 * pkin(6) + t9;
t3 = qJD(2) * pkin(3) - t16 * pkin(6) + t8;
t2 = t27 * t3 + t40 * t4;
t1 = -t27 * t4 + t40 * t3;
t6 = [0, 0, 0, 0, 0, t41, 0, 0, 0, 0, t24 * t41, t28 * t39, t33, t25 * t41, t32, t23, pkin(1) * t39 - pkin(5) * t33, -t30 * pkin(1) * t28 - pkin(5) * t32, (t24 + t25) * t30 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t25 / 0.2e1 + t24 / 0.2e1) * pkin(5) ^ 2) * t30, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * qJD(2), t14 ^ 2 / 0.2e1, -t14 * qJD(2), t23, t8 * qJD(2) + t20 * t14, -t9 * qJD(2) + t20 * t16, -t9 * t14 - t8 * t16, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t7 ^ 2 / 0.2e1, -t7 * t5, t7 * t22, t5 ^ 2 / 0.2e1, -t5 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t10 * t5, t10 * t7 - t2 * t22, -t1 * t7 - t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t6;
