% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP3
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:12
% EndTime: 2019-12-31 19:51:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (274->34), mult. (430->84), div. (0->0), fcn. (233->6), ass. (0->39)
t25 = sin(pkin(8));
t21 = t25 ^ 2;
t44 = t21 / 0.2e1;
t26 = cos(pkin(8));
t22 = t26 ^ 2;
t43 = t22 / 0.2e1;
t27 = sin(qJ(4));
t42 = cos(qJ(4));
t23 = qJD(1) + qJD(2);
t28 = sin(qJ(2));
t38 = pkin(1) * qJD(1);
t35 = t28 * t38;
t18 = t23 * qJ(3) + t35;
t33 = pkin(7) * t23 + t18;
t8 = t33 * t25;
t9 = t33 * t26;
t5 = -t27 * t8 + t42 * t9;
t39 = t23 * t26;
t40 = t23 * t25;
t13 = t27 * t40 - t42 * t39;
t15 = (t42 * t25 + t26 * t27) * t23;
t41 = t15 * t13;
t37 = qJD(4) * t13;
t36 = t13 ^ 2 / 0.2e1;
t29 = cos(qJ(2));
t34 = t29 * t38;
t32 = qJD(3) - t34;
t4 = -t27 * t9 - t42 * t8;
t12 = (-pkin(3) * t26 - pkin(2)) * t23 + t32;
t30 = qJD(1) ^ 2;
t24 = qJD(4) ^ 2 / 0.2e1;
t20 = t23 ^ 2;
t16 = -t23 * pkin(2) + t32;
t11 = t15 * qJD(4);
t10 = t15 ^ 2 / 0.2e1;
t3 = qJD(4) * qJ(5) + t5;
t2 = -qJD(4) * pkin(4) + qJD(5) - t4;
t1 = t13 * pkin(4) - t15 * qJ(5) + t12;
t6 = [0, 0, 0, 0, 0, t30 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 / 0.2e1, t23 * t34, -t23 * t35, 0, (t28 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t30, t20 * t44, t25 * t20 * t26, 0, t20 * t43, 0, 0, -t16 * t39, t16 * t40, (t21 + t22) * t23 * t18, t16 ^ 2 / 0.2e1 + (t43 + t44) * t18 ^ 2, t10, -t41, t11, t36, -t37, t24, t4 * qJD(4) + t12 * t13, -t5 * qJD(4) + t12 * t15, -t5 * t13 - t4 * t15, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10, t11, t41, t24, t37, t36, -t2 * qJD(4) + t1 * t13, -t3 * t13 + t2 * t15, t3 * qJD(4) - t1 * t15, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t6;
