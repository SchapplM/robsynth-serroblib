% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:17:46
% EndTime: 2019-03-09 00:17:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (511->57), mult. (1149->132), div. (0->0), fcn. (815->10), ass. (0->50)
t50 = qJD(2) ^ 2;
t64 = t50 / 0.2e1;
t44 = sin(qJ(5));
t62 = cos(qJ(5));
t47 = sin(qJ(2));
t42 = sin(pkin(6));
t59 = qJD(1) * t42;
t32 = qJD(2) * pkin(8) + t47 * t59;
t46 = sin(qJ(3));
t48 = cos(qJ(3));
t43 = cos(pkin(6));
t58 = qJD(1) * t43;
t24 = t48 * t32 + t46 * t58;
t20 = qJD(3) * pkin(9) + t24;
t49 = cos(qJ(2));
t53 = t49 * t59;
t25 = -t53 + (-pkin(3) * t48 - pkin(9) * t46 - pkin(2)) * qJD(2);
t45 = sin(qJ(4));
t63 = cos(qJ(4));
t10 = -t45 * t20 + t63 * t25;
t57 = qJD(2) * t46;
t31 = t45 * qJD(3) + t63 * t57;
t56 = t48 * qJD(2);
t38 = -qJD(4) + t56;
t7 = -t38 * pkin(4) - t31 * pkin(10) + t10;
t11 = t63 * t20 + t45 * t25;
t29 = -t63 * qJD(3) + t45 * t57;
t9 = -t29 * pkin(10) + t11;
t4 = t44 * t7 + t62 * t9;
t15 = t62 * t29 + t44 * t31;
t17 = -t44 * t29 + t62 * t31;
t61 = t17 * t15;
t35 = -qJD(5) + t38;
t60 = t35 * t15;
t55 = t15 ^ 2 / 0.2e1;
t54 = qJD(2) * qJD(3);
t52 = qJD(2) * t59;
t23 = -t46 * t32 + t48 * t58;
t3 = -t44 * t9 + t62 * t7;
t19 = -qJD(3) * pkin(3) - t23;
t13 = t29 * pkin(4) + t19;
t51 = qJD(1) ^ 2;
t34 = t35 ^ 2 / 0.2e1;
t33 = -qJD(2) * pkin(2) - t53;
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t35;
t5 = t15 * pkin(5) - t17 * qJ(6) + t13;
t2 = -t35 * qJ(6) + t4;
t1 = t35 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t51 / 0.2e1, 0, 0, 0, 0, 0, t64, t49 * t52, -t47 * t52, 0 (t43 ^ 2 / 0.2e1 + (t47 ^ 2 / 0.2e1 + t49 ^ 2 / 0.2e1) * t42 ^ 2) * t51, t46 ^ 2 * t64, t46 * t50 * t48, t46 * t54, t48 ^ 2 * t64, t48 * t54, qJD(3) ^ 2 / 0.2e1, t23 * qJD(3) - t33 * t56, -t24 * qJD(3) + t33 * t57 (-t23 * t46 + t24 * t48) * qJD(2), t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t38, t29 ^ 2 / 0.2e1, t29 * t38, t38 ^ 2 / 0.2e1, -t10 * t38 + t19 * t29, t11 * t38 + t19 * t31, -t10 * t31 - t11 * t29, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t14, -t61, -t12, t55, t60, t34, t13 * t15 - t3 * t35, t13 * t17 + t4 * t35, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t14, -t12, t61, t34, -t60, t55, t1 * t35 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 - t2 * t35, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
