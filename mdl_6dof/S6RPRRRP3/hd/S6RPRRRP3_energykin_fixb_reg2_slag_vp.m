% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:15
% EndTime: 2019-03-09 06:05:15
% DurationCPUTime: 0.16s
% Computational Cost: add. (511->56), mult. (1102->128), div. (0->0), fcn. (687->8), ass. (0->45)
t48 = qJD(1) ^ 2;
t41 = t48 / 0.2e1;
t44 = sin(qJ(5));
t57 = cos(qJ(5));
t42 = sin(pkin(10));
t33 = (pkin(1) * t42 + pkin(7)) * qJD(1);
t46 = sin(qJ(3));
t47 = cos(qJ(3));
t25 = t46 * qJD(2) + t47 * t33;
t22 = qJD(3) * pkin(8) + t25;
t43 = cos(pkin(10));
t50 = -pkin(1) * t43 - pkin(2);
t23 = (-pkin(3) * t47 - pkin(8) * t46 + t50) * qJD(1);
t45 = sin(qJ(4));
t58 = cos(qJ(4));
t10 = -t45 * t22 + t58 * t23;
t54 = qJD(1) * t46;
t31 = t45 * qJD(3) + t58 * t54;
t53 = t47 * qJD(1);
t37 = -qJD(4) + t53;
t7 = -t37 * pkin(4) - t31 * pkin(9) + t10;
t11 = t58 * t22 + t45 * t23;
t29 = -t58 * qJD(3) + t45 * t54;
t9 = -t29 * pkin(9) + t11;
t4 = t44 * t7 + t57 * t9;
t59 = pkin(1) * t48;
t15 = t57 * t29 + t44 * t31;
t17 = -t44 * t29 + t57 * t31;
t56 = t17 * t15;
t35 = -qJD(5) + t37;
t55 = t35 * t15;
t52 = t15 ^ 2 / 0.2e1;
t51 = qJD(1) * qJD(3);
t24 = t47 * qJD(2) - t46 * t33;
t3 = -t44 * t9 + t57 * t7;
t21 = -qJD(3) * pkin(3) - t24;
t13 = t29 * pkin(4) + t21;
t34 = t50 * qJD(1);
t32 = t35 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t35;
t5 = t15 * pkin(5) - t17 * qJ(6) + t13;
t2 = -t35 * qJ(6) + t4;
t1 = t35 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t43 * t59, -t42 * t59, 0, qJD(2) ^ 2 / 0.2e1 + (t42 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t48, t46 ^ 2 * t41, t47 * t48 * t46, t46 * t51, t47 ^ 2 * t41, t47 * t51, qJD(3) ^ 2 / 0.2e1, t24 * qJD(3) - t34 * t53, -t25 * qJD(3) + t34 * t54 (-t24 * t46 + t25 * t47) * qJD(1), t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t37, t29 ^ 2 / 0.2e1, t29 * t37, t37 ^ 2 / 0.2e1, -t10 * t37 + t21 * t29, t11 * t37 + t21 * t31, -t10 * t31 - t11 * t29, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14, -t56, -t12, t52, t55, t32, t13 * t15 - t3 * t35, t13 * t17 + t4 * t35, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t14, -t12, t56, t32, -t55, t52, t1 * t35 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 - t2 * t35, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
