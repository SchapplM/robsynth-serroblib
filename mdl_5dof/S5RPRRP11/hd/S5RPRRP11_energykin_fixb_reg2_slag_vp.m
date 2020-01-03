% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:31
% EndTime: 2019-12-31 18:54:31
% DurationCPUTime: 0.17s
% Computational Cost: add. (311->45), mult. (824->102), div. (0->0), fcn. (549->6), ass. (0->39)
t38 = qJD(1) ^ 2;
t48 = t38 / 0.2e1;
t34 = sin(pkin(8));
t42 = qJD(1) * t34;
t43 = pkin(6) + qJ(2);
t26 = t43 * t42;
t35 = cos(pkin(8));
t41 = qJD(1) * t35;
t27 = t43 * t41;
t37 = sin(qJ(3));
t47 = cos(qJ(3));
t13 = -t37 * t26 + t47 * t27;
t11 = qJD(3) * pkin(7) + t13;
t36 = sin(qJ(4));
t46 = cos(qJ(4));
t23 = t37 * t42 - t47 * t41;
t25 = (t47 * t34 + t35 * t37) * qJD(1);
t28 = qJD(2) + (-pkin(2) * t35 - pkin(1)) * qJD(1);
t7 = t23 * pkin(3) - t25 * pkin(7) + t28;
t4 = t46 * t11 + t36 * t7;
t15 = -t46 * qJD(3) + t36 * t25;
t17 = t36 * qJD(3) + t46 * t25;
t45 = t17 * t15;
t19 = qJD(4) + t23;
t44 = t19 * t15;
t40 = t15 ^ 2 / 0.2e1;
t12 = -t47 * t26 - t37 * t27;
t3 = -t36 * t11 + t46 * t7;
t10 = -qJD(3) * pkin(3) - t12;
t33 = t35 ^ 2;
t32 = t34 ^ 2;
t30 = -qJD(1) * pkin(1) + qJD(2);
t18 = t19 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t8 = t17 * t19;
t5 = t15 * pkin(4) - t17 * qJ(5) + t10;
t2 = t19 * qJ(5) + t4;
t1 = -t19 * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t48, 0, 0, 0, 0, t32 * t48, t34 * t38 * t35, 0, t33 * t48, 0, 0, -t30 * t41, t30 * t42, (t32 + t33) * t38 * qJ(2), t30 ^ 2 / 0.2e1 + (t33 / 0.2e1 + t32 / 0.2e1) * qJ(2) ^ 2 * t38, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * qJD(3), t23 ^ 2 / 0.2e1, -t23 * qJD(3), qJD(3) ^ 2 / 0.2e1, t12 * qJD(3) + t28 * t23, -t13 * qJD(3) + t28 * t25, -t12 * t25 - t13 * t23, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t14, -t45, t8, t40, -t44, t18, t10 * t15 + t3 * t19, t10 * t17 - t4 * t19, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t14, t8, t45, t18, t44, t40, -t1 * t19 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t19, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
