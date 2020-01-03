% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP9
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:27
% EndTime: 2019-12-31 20:07:27
% DurationCPUTime: 0.15s
% Computational Cost: add. (363->49), mult. (885->109), div. (0->0), fcn. (558->6), ass. (0->42)
t40 = qJD(1) ^ 2;
t53 = t40 / 0.2e1;
t37 = sin(qJ(4));
t52 = cos(qJ(4));
t38 = sin(qJ(2));
t39 = cos(qJ(2));
t21 = (-pkin(2) * t39 - qJ(3) * t38 - pkin(1)) * qJD(1);
t46 = t39 * qJD(1);
t27 = pkin(6) * t46 + qJD(2) * qJ(3);
t36 = sin(pkin(8));
t48 = cos(pkin(8));
t15 = t48 * t21 - t36 * t27;
t47 = qJD(1) * t38;
t24 = t36 * qJD(2) + t48 * t47;
t7 = -pkin(3) * t46 - t24 * pkin(7) + t15;
t16 = t36 * t21 + t48 * t27;
t22 = -t48 * qJD(2) + t36 * t47;
t9 = -t22 * pkin(7) + t16;
t4 = t37 * t7 + t52 * t9;
t12 = t52 * t22 + t37 * t24;
t14 = -t37 * t22 + t52 * t24;
t51 = t14 * t12;
t29 = -qJD(4) + t46;
t50 = t29 * t12;
t49 = t39 * t40;
t45 = t12 ^ 2 / 0.2e1;
t44 = qJD(1) * qJD(2);
t43 = t38 * t44;
t42 = t39 * t44;
t26 = -qJD(2) * pkin(2) + pkin(6) * t47 + qJD(3);
t3 = -t37 * t9 + t52 * t7;
t17 = t22 * pkin(3) + t26;
t35 = t39 ^ 2;
t34 = t38 ^ 2;
t31 = t35 * t53;
t28 = t29 ^ 2 / 0.2e1;
t11 = t14 ^ 2 / 0.2e1;
t10 = t14 * t29;
t5 = t12 * pkin(4) - t14 * qJ(5) + t17;
t2 = -t29 * qJ(5) + t4;
t1 = t29 * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t53, 0, 0, 0, 0, t34 * t53, t38 * t49, t43, t31, t42, qJD(2) ^ 2 / 0.2e1, pkin(1) * t49 - pkin(6) * t43, -t40 * pkin(1) * t38 - pkin(6) * t42, (t34 + t35) * t40 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t35 / 0.2e1 + t34 / 0.2e1) * pkin(6) ^ 2) * t40, t24 ^ 2 / 0.2e1, -t24 * t22, -t24 * t46, t22 ^ 2 / 0.2e1, t22 * t46, t31, -t15 * t46 + t26 * t22, t16 * t46 + t26 * t24, -t15 * t24 - t16 * t22, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t11, -t51, -t10, t45, t50, t28, t17 * t12 - t3 * t29, t17 * t14 + t4 * t29, -t4 * t12 - t3 * t14, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t11, -t10, t51, t28, -t50, t45, t1 * t29 + t5 * t12, t1 * t14 - t2 * t12, -t5 * t14 - t2 * t29, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
