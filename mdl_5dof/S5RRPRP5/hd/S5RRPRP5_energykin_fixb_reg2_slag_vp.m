% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP5
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:55
% EndTime: 2019-12-31 19:54:56
% DurationCPUTime: 0.17s
% Computational Cost: add. (361->46), mult. (973->107), div. (0->0), fcn. (647->6), ass. (0->43)
t38 = qJD(1) ^ 2;
t52 = t38 / 0.2e1;
t35 = sin(qJ(4));
t51 = cos(qJ(4));
t36 = sin(qJ(2));
t45 = qJD(1) * t36;
t47 = pkin(6) + qJ(3);
t25 = qJD(2) * pkin(2) - t47 * t45;
t37 = cos(qJ(2));
t44 = qJD(1) * t37;
t26 = t47 * t44;
t34 = sin(pkin(8));
t46 = cos(pkin(8));
t15 = t46 * t25 - t34 * t26;
t23 = (t34 * t37 + t46 * t36) * qJD(1);
t8 = qJD(2) * pkin(3) - t23 * pkin(7) + t15;
t16 = t34 * t25 + t46 * t26;
t21 = t34 * t45 - t46 * t44;
t9 = -t21 * pkin(7) + t16;
t4 = t35 * t8 + t51 * t9;
t12 = t51 * t21 + t35 * t23;
t14 = -t35 * t21 + t51 * t23;
t50 = t14 * t12;
t30 = qJD(2) + qJD(4);
t49 = t30 * t12;
t48 = t37 * t38;
t43 = t12 ^ 2 / 0.2e1;
t42 = qJD(1) * qJD(2);
t41 = t36 * t42;
t40 = t37 * t42;
t3 = -t35 * t9 + t51 * t8;
t27 = qJD(3) + (-pkin(2) * t37 - pkin(1)) * qJD(1);
t17 = t21 * pkin(3) + t27;
t33 = t37 ^ 2;
t32 = t36 ^ 2;
t31 = qJD(2) ^ 2 / 0.2e1;
t29 = t30 ^ 2 / 0.2e1;
t11 = t14 ^ 2 / 0.2e1;
t10 = t14 * t30;
t5 = t12 * pkin(4) - t14 * qJ(5) + t17;
t2 = t30 * qJ(5) + t4;
t1 = -t30 * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t52, 0, 0, 0, 0, t32 * t52, t36 * t48, t41, t33 * t52, t40, t31, pkin(1) * t48 - pkin(6) * t41, -t38 * pkin(1) * t36 - pkin(6) * t40, (t32 + t33) * t38 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t33 / 0.2e1 + t32 / 0.2e1) * pkin(6) ^ 2) * t38, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * qJD(2), t21 ^ 2 / 0.2e1, -t21 * qJD(2), t31, t15 * qJD(2) + t27 * t21, -t16 * qJD(2) + t27 * t23, -t15 * t23 - t16 * t21, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t11, -t50, t10, t43, -t49, t29, t17 * t12 + t3 * t30, t17 * t14 - t4 * t30, -t4 * t12 - t3 * t14, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t11, t10, t50, t29, t49, t43, -t1 * t30 + t5 * t12, t1 * t14 - t2 * t12, -t5 * t14 + t2 * t30, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
