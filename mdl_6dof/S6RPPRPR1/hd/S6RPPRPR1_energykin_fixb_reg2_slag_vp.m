% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:01
% EndTime: 2019-03-09 01:40:01
% DurationCPUTime: 0.16s
% Computational Cost: add. (587->60), mult. (1393->139), div. (0->0), fcn. (955->10), ass. (0->45)
t50 = qJD(1) ^ 2;
t42 = t50 / 0.2e1;
t59 = pkin(1) * t50;
t45 = sin(pkin(9));
t37 = (pkin(1) * t45 + qJ(3)) * qJD(1);
t46 = cos(pkin(10));
t41 = t46 * qJD(2);
t44 = sin(pkin(10));
t26 = t41 + (-pkin(7) * qJD(1) - t37) * t44;
t29 = t44 * qJD(2) + t46 * t37;
t54 = qJD(1) * t46;
t27 = pkin(7) * t54 + t29;
t49 = sin(qJ(4));
t58 = cos(qJ(4));
t15 = t49 * t26 + t58 * t27;
t10 = qJD(4) * qJ(5) + t15;
t47 = cos(pkin(9));
t52 = -pkin(1) * t47 - pkin(2);
t31 = qJD(3) + (-pkin(3) * t46 + t52) * qJD(1);
t55 = qJD(1) * t44;
t32 = t49 * t55 - t58 * t54;
t34 = (t58 * t44 + t46 * t49) * qJD(1);
t18 = t32 * pkin(4) - t34 * qJ(5) + t31;
t43 = sin(pkin(11));
t56 = cos(pkin(11));
t6 = t56 * t10 + t43 * t18;
t57 = cos(qJ(6));
t53 = t32 ^ 2 / 0.2e1;
t5 = -t43 * t10 + t56 * t18;
t14 = t58 * t26 - t49 * t27;
t9 = -qJD(4) * pkin(4) + qJD(5) - t14;
t48 = sin(qJ(6));
t36 = t52 * qJD(1) + qJD(3);
t30 = qJD(6) + t32;
t28 = -t44 * t37 + t41;
t25 = t43 * qJD(4) + t56 * t34;
t23 = -t56 * qJD(4) + t43 * t34;
t13 = -t48 * t23 + t57 * t25;
t11 = t57 * t23 + t48 * t25;
t7 = t23 * pkin(5) + t9;
t4 = -t23 * pkin(8) + t6;
t3 = t32 * pkin(5) - t25 * pkin(8) + t5;
t2 = t48 * t3 + t57 * t4;
t1 = t57 * t3 - t48 * t4;
t8 = [0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t47 * t59, -t45 * t59, 0, qJD(2) ^ 2 / 0.2e1 + (t45 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t50, t44 ^ 2 * t42, t44 * t50 * t46, 0, t46 ^ 2 * t42, 0, 0, -t36 * t54, t36 * t55 (-t28 * t44 + t29 * t46) * qJD(1), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * qJD(4), t53, -t32 * qJD(4), qJD(4) ^ 2 / 0.2e1, t14 * qJD(4) + t31 * t32, -t15 * qJD(4) + t31 * t34, -t14 * t34 - t15 * t32, t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t32, t23 ^ 2 / 0.2e1, -t23 * t32, t53, t9 * t23 + t5 * t32, t9 * t25 - t6 * t32, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t30, t11 ^ 2 / 0.2e1, -t11 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t11, t7 * t13 - t2 * t30, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
