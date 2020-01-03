% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:10
% EndTime: 2019-12-31 18:49:10
% DurationCPUTime: 0.17s
% Computational Cost: add. (333->45), mult. (930->102), div. (0->0), fcn. (643->6), ass. (0->39)
t38 = qJD(1) ^ 2;
t48 = t38 / 0.2e1;
t36 = sin(qJ(4));
t46 = cos(qJ(4));
t34 = sin(pkin(8));
t42 = qJD(1) * t34;
t43 = pkin(6) + qJ(2);
t25 = t43 * t42;
t35 = cos(pkin(8));
t41 = qJD(1) * t35;
t26 = t43 * t41;
t37 = sin(qJ(3));
t47 = cos(qJ(3));
t15 = -t47 * t25 - t26 * t37;
t24 = (t34 * t47 + t35 * t37) * qJD(1);
t8 = qJD(3) * pkin(3) - pkin(7) * t24 + t15;
t16 = -t37 * t25 + t47 * t26;
t22 = t37 * t42 - t41 * t47;
t9 = -pkin(7) * t22 + t16;
t4 = t36 * t8 + t46 * t9;
t12 = t22 * t46 + t36 * t24;
t14 = -t36 * t22 + t24 * t46;
t45 = t14 * t12;
t33 = qJD(3) + qJD(4);
t44 = t33 * t12;
t40 = t12 ^ 2 / 0.2e1;
t3 = -t36 * t9 + t46 * t8;
t27 = qJD(2) + (-pkin(2) * t35 - pkin(1)) * qJD(1);
t17 = t22 * pkin(3) + t27;
t32 = t35 ^ 2;
t31 = t34 ^ 2;
t30 = -qJD(1) * pkin(1) + qJD(2);
t29 = t33 ^ 2 / 0.2e1;
t11 = t14 ^ 2 / 0.2e1;
t10 = t14 * t33;
t5 = t12 * pkin(4) - t14 * qJ(5) + t17;
t2 = qJ(5) * t33 + t4;
t1 = -t33 * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t48, 0, 0, 0, 0, t31 * t48, t34 * t38 * t35, 0, t32 * t48, 0, 0, -t30 * t41, t30 * t42, (t31 + t32) * t38 * qJ(2), t30 ^ 2 / 0.2e1 + (t32 / 0.2e1 + t31 / 0.2e1) * qJ(2) ^ 2 * t38, t24 ^ 2 / 0.2e1, -t24 * t22, t24 * qJD(3), t22 ^ 2 / 0.2e1, -t22 * qJD(3), qJD(3) ^ 2 / 0.2e1, qJD(3) * t15 + t22 * t27, -qJD(3) * t16 + t24 * t27, -t15 * t24 - t16 * t22, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t11, -t45, t10, t40, -t44, t29, t12 * t17 + t3 * t33, t14 * t17 - t33 * t4, -t12 * t4 - t14 * t3, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t11, t10, t45, t29, t44, t40, -t1 * t33 + t12 * t5, t1 * t14 - t12 * t2, -t14 * t5 + t2 * t33, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
