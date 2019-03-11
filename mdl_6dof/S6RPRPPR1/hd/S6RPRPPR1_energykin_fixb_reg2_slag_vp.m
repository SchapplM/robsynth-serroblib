% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:33
% EndTime: 2019-03-09 02:39:33
% DurationCPUTime: 0.16s
% Computational Cost: add. (619->61), mult. (1441->144), div. (0->0), fcn. (963->10), ass. (0->48)
t51 = qJD(1) ^ 2;
t43 = t51 / 0.2e1;
t62 = pkin(1) * t51;
t46 = sin(pkin(9));
t36 = (pkin(1) * t46 + pkin(7)) * qJD(1);
t50 = cos(qJ(3));
t41 = t50 * qJD(2);
t49 = sin(qJ(3));
t55 = qJ(4) * qJD(1);
t23 = qJD(3) * pkin(3) + t41 + (-t36 - t55) * t49;
t29 = t49 * qJD(2) + t50 * t36;
t27 = t50 * t55 + t29;
t45 = sin(pkin(10));
t60 = cos(pkin(10));
t12 = t45 * t23 + t60 * t27;
t10 = qJD(3) * qJ(5) + t12;
t47 = cos(pkin(9));
t53 = -pkin(1) * t47 - pkin(2);
t31 = qJD(4) + (-pkin(3) * t50 + t53) * qJD(1);
t57 = qJD(1) * t50;
t58 = qJD(1) * t49;
t32 = t45 * t58 - t60 * t57;
t34 = (t45 * t50 + t60 * t49) * qJD(1);
t18 = t32 * pkin(4) - t34 * qJ(5) + t31;
t44 = sin(pkin(11));
t59 = cos(pkin(11));
t6 = t59 * t10 + t44 * t18;
t61 = cos(qJ(6));
t56 = t32 ^ 2 / 0.2e1;
t54 = qJD(1) * qJD(3);
t5 = -t44 * t10 + t59 * t18;
t11 = t60 * t23 - t45 * t27;
t9 = -qJD(3) * pkin(4) + qJD(5) - t11;
t48 = sin(qJ(6));
t42 = qJD(3) ^ 2 / 0.2e1;
t37 = t53 * qJD(1);
t30 = qJD(6) + t32;
t28 = -t49 * t36 + t41;
t26 = t44 * qJD(3) + t59 * t34;
t24 = -t59 * qJD(3) + t44 * t34;
t15 = -t48 * t24 + t61 * t26;
t13 = t61 * t24 + t48 * t26;
t7 = t24 * pkin(5) + t9;
t4 = -t24 * pkin(8) + t6;
t3 = t32 * pkin(5) - t26 * pkin(8) + t5;
t2 = t48 * t3 + t61 * t4;
t1 = t61 * t3 - t48 * t4;
t8 = [0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t47 * t62, -t46 * t62, 0, qJD(2) ^ 2 / 0.2e1 + (t46 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t51, t49 ^ 2 * t43, t49 * t51 * t50, t49 * t54, t50 ^ 2 * t43, t50 * t54, t42, t28 * qJD(3) - t37 * t57, -t29 * qJD(3) + t37 * t58 (-t28 * t49 + t29 * t50) * qJD(1), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t32 * t34, qJD(3) * t34, t56, -t32 * qJD(3), t42, t11 * qJD(3) + t31 * t32, -t12 * qJD(3) + t31 * t34, -t11 * t34 - t12 * t32, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, t26 * t32, t24 ^ 2 / 0.2e1, -t24 * t32, t56, t9 * t24 + t5 * t32, t9 * t26 - t6 * t32, -t6 * t24 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t30, t13 ^ 2 / 0.2e1, -t13 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t13, t7 * t15 - t2 * t30, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
