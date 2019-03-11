% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:21
% EndTime: 2019-03-09 02:21:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (597->60), mult. (1393->142), div. (0->0), fcn. (955->10), ass. (0->45)
t52 = qJD(1) ^ 2;
t43 = t52 / 0.2e1;
t59 = pkin(1) * t52;
t45 = sin(pkin(10));
t38 = (pkin(1) * t45 + qJ(3)) * qJD(1);
t46 = cos(pkin(11));
t41 = t46 * qJD(2);
t44 = sin(pkin(11));
t23 = t41 + (-pkin(7) * qJD(1) - t38) * t44;
t29 = t44 * qJD(2) + t46 * t38;
t55 = qJD(1) * t46;
t24 = pkin(7) * t55 + t29;
t50 = sin(qJ(4));
t51 = cos(qJ(4));
t12 = t50 * t23 + t51 * t24;
t10 = qJD(4) * pkin(8) + t12;
t47 = cos(pkin(10));
t54 = -pkin(1) * t47 - pkin(2);
t32 = qJD(3) + (-pkin(3) * t46 + t54) * qJD(1);
t56 = qJD(1) * t44;
t33 = t50 * t56 - t51 * t55;
t35 = (t44 * t51 + t46 * t50) * qJD(1);
t18 = t33 * pkin(4) - t35 * pkin(8) + t32;
t49 = sin(qJ(5));
t58 = cos(qJ(5));
t6 = t58 * t10 + t49 * t18;
t57 = cos(qJ(6));
t5 = -t49 * t10 + t58 * t18;
t11 = t51 * t23 - t50 * t24;
t31 = qJD(5) + t33;
t9 = -qJD(4) * pkin(4) - t11;
t48 = sin(qJ(6));
t37 = t54 * qJD(1) + qJD(3);
t30 = qJD(6) + t31;
t28 = -t44 * t38 + t41;
t27 = t49 * qJD(4) + t58 * t35;
t25 = -t58 * qJD(4) + t49 * t35;
t15 = -t48 * t25 + t57 * t27;
t13 = t57 * t25 + t48 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(9) + t6;
t3 = t31 * pkin(5) - t27 * pkin(9) + t5;
t2 = t48 * t3 + t57 * t4;
t1 = t57 * t3 - t48 * t4;
t8 = [0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t47 * t59, -t45 * t59, 0, qJD(2) ^ 2 / 0.2e1 + (t45 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t52, t44 ^ 2 * t43, t44 * t52 * t46, 0, t46 ^ 2 * t43, 0, 0, -t37 * t55, t37 * t56 (-t28 * t44 + t29 * t46) * qJD(1), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t33 * t35, qJD(4) * t35, t33 ^ 2 / 0.2e1, -t33 * qJD(4), qJD(4) ^ 2 / 0.2e1, t11 * qJD(4) + t32 * t33, -t12 * qJD(4) + t32 * t35, -t11 * t35 - t12 * t33, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t31, t25 ^ 2 / 0.2e1, -t25 * t31, t31 ^ 2 / 0.2e1, t9 * t25 + t5 * t31, t9 * t27 - t6 * t31, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t30, t13 ^ 2 / 0.2e1, -t13 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t13, t7 * t15 - t2 * t30, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
