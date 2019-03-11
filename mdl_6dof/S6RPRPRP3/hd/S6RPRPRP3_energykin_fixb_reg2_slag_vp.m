% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:35
% EndTime: 2019-03-09 03:09:36
% DurationCPUTime: 0.14s
% Computational Cost: add. (489->56), mult. (1103->125), div. (0->0), fcn. (687->8), ass. (0->45)
t47 = qJD(1) ^ 2;
t40 = t47 / 0.2e1;
t44 = sin(qJ(5));
t57 = cos(qJ(5));
t42 = sin(pkin(9));
t32 = (pkin(1) * t42 + pkin(7)) * qJD(1);
t45 = sin(qJ(3));
t46 = cos(qJ(3));
t25 = t45 * qJD(2) + t46 * t32;
t22 = qJD(3) * qJ(4) + t25;
t43 = cos(pkin(9));
t49 = -pkin(1) * t43 - pkin(2);
t23 = (-pkin(3) * t46 - qJ(4) * t45 + t49) * qJD(1);
t41 = sin(pkin(10));
t54 = cos(pkin(10));
t10 = -t41 * t22 + t54 * t23;
t53 = qJD(1) * t45;
t29 = t41 * qJD(3) + t54 * t53;
t52 = t46 * qJD(1);
t7 = -pkin(4) * t52 - t29 * pkin(8) + t10;
t11 = t54 * t22 + t41 * t23;
t27 = -t54 * qJD(3) + t41 * t53;
t9 = -t27 * pkin(8) + t11;
t4 = t44 * t7 + t57 * t9;
t58 = pkin(1) * t47;
t15 = t57 * t27 + t44 * t29;
t17 = -t44 * t27 + t57 * t29;
t56 = t17 * t15;
t35 = -qJD(5) + t52;
t55 = t35 * t15;
t51 = t15 ^ 2 / 0.2e1;
t50 = qJD(1) * qJD(3);
t24 = t46 * qJD(2) - t45 * t32;
t3 = -t44 * t9 + t57 * t7;
t21 = -qJD(3) * pkin(3) + qJD(4) - t24;
t13 = t27 * pkin(4) + t21;
t37 = t46 ^ 2 * t40;
t34 = t35 ^ 2 / 0.2e1;
t33 = t49 * qJD(1);
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t35;
t5 = t15 * pkin(5) - t17 * qJ(6) + t13;
t2 = -t35 * qJ(6) + t4;
t1 = t35 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43 * t58, -t42 * t58, 0, qJD(2) ^ 2 / 0.2e1 + (t42 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t47, t45 ^ 2 * t40, t45 * t47 * t46, t45 * t50, t37, t46 * t50, qJD(3) ^ 2 / 0.2e1, t24 * qJD(3) - t33 * t52, -t25 * qJD(3) + t33 * t53 (-t24 * t45 + t25 * t46) * qJD(1), t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, -t29 * t52, t27 ^ 2 / 0.2e1, t27 * t52, t37, -t10 * t52 + t21 * t27, t11 * t52 + t21 * t29, -t10 * t29 - t11 * t27, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14, -t56, -t12, t51, t55, t34, t13 * t15 - t3 * t35, t13 * t17 + t4 * t35, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t14, -t12, t56, t34, -t55, t51, t1 * t35 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 - t2 * t35, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
