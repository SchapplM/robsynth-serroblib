% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP1
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:05
% EndTime: 2019-03-09 05:57:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (496->55), mult. (1076->130), div. (0->0), fcn. (685->8), ass. (0->46)
t46 = qJD(1) ^ 2;
t39 = t46 / 0.2e1;
t57 = pkin(1) * t46;
t43 = sin(qJ(4));
t45 = cos(qJ(3));
t51 = qJD(1) * t45;
t44 = sin(qJ(3));
t52 = qJD(1) * t44;
t56 = cos(qJ(4));
t27 = t43 * t52 - t56 * t51;
t29 = (t43 * t45 + t56 * t44) * qJD(1);
t41 = cos(pkin(10));
t48 = -pkin(1) * t41 - pkin(2);
t30 = (-pkin(3) * t45 + t48) * qJD(1);
t12 = t27 * pkin(4) - t29 * pkin(9) + t30;
t42 = sin(qJ(5));
t55 = cos(qJ(5));
t40 = sin(pkin(10));
t32 = (pkin(1) * t40 + pkin(7)) * qJD(1);
t37 = t45 * qJD(2);
t18 = qJD(3) * pkin(3) + t37 + (-pkin(8) * qJD(1) - t32) * t44;
t25 = t44 * qJD(2) + t45 * t32;
t22 = pkin(8) * t51 + t25;
t10 = t43 * t18 + t56 * t22;
t38 = qJD(3) + qJD(4);
t8 = t38 * pkin(9) + t10;
t4 = t42 * t12 + t55 * t8;
t19 = t42 * t29 - t55 * t38;
t21 = t55 * t29 + t42 * t38;
t54 = t21 * t19;
t26 = qJD(5) + t27;
t53 = t26 * t19;
t50 = t19 ^ 2 / 0.2e1;
t49 = qJD(1) * qJD(3);
t9 = t56 * t18 - t43 * t22;
t3 = t55 * t12 - t42 * t8;
t7 = -t38 * pkin(4) - t9;
t33 = t48 * qJD(1);
t24 = -t44 * t32 + t37;
t23 = t26 ^ 2 / 0.2e1;
t15 = t21 ^ 2 / 0.2e1;
t13 = t21 * t26;
t5 = t19 * pkin(5) - t21 * qJ(6) + t7;
t2 = t26 * qJ(6) + t4;
t1 = -t26 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t41 * t57, -t40 * t57, 0, qJD(2) ^ 2 / 0.2e1 + (t40 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t46, t44 ^ 2 * t39, t44 * t46 * t45, t44 * t49, t45 ^ 2 * t39, t45 * t49, qJD(3) ^ 2 / 0.2e1, t24 * qJD(3) - t33 * t51, -t25 * qJD(3) + t33 * t52 (-t24 * t44 + t25 * t45) * qJD(1), t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t38, t27 ^ 2 / 0.2e1, -t27 * t38, t38 ^ 2 / 0.2e1, t30 * t27 + t9 * t38, -t10 * t38 + t30 * t29, -t10 * t27 - t9 * t29, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t15, -t54, t13, t50, -t53, t23, t7 * t19 + t3 * t26, t7 * t21 - t4 * t26, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t15, t13, t54, t23, t53, t50, -t1 * t26 + t5 * t19, t1 * t21 - t2 * t19, t2 * t26 - t5 * t21, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
