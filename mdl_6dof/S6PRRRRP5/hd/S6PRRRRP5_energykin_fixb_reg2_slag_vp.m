% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:25:52
% EndTime: 2019-03-09 00:25:52
% DurationCPUTime: 0.23s
% Computational Cost: add. (727->60), mult. (1776->140), div. (0->0), fcn. (1412->12), ass. (0->56)
t59 = cos(qJ(2));
t51 = sin(pkin(6));
t70 = qJD(1) * t51;
t41 = qJD(2) * pkin(2) + t59 * t70;
t50 = sin(pkin(7));
t52 = cos(pkin(7));
t53 = cos(pkin(6));
t69 = qJD(1) * t53;
t75 = t41 * t52 + t50 * t69;
t57 = sin(qJ(2));
t68 = qJD(2) * t50;
t39 = pkin(9) * t68 + t57 * t70;
t56 = sin(qJ(3));
t58 = cos(qJ(3));
t21 = -t39 * t56 + t58 * t75;
t47 = qJD(2) * t52 + qJD(3);
t19 = -t47 * pkin(3) - t21;
t55 = sin(qJ(4));
t66 = t56 * t68;
t74 = cos(qJ(4));
t33 = -t47 * t74 + t55 * t66;
t35 = t47 * t55 + t66 * t74;
t13 = t33 * pkin(4) - t35 * pkin(11) + t19;
t54 = sin(qJ(5));
t73 = cos(qJ(5));
t22 = t58 * t39 + t56 * t75;
t20 = pkin(10) * t47 + t22;
t46 = t52 * t69;
t24 = t46 + (-t41 + (-pkin(3) * t58 - pkin(10) * t56) * qJD(2)) * t50;
t12 = t20 * t74 + t24 * t55;
t65 = t58 * t68;
t44 = -qJD(4) + t65;
t8 = -pkin(11) * t44 + t12;
t4 = t13 * t54 + t73 * t8;
t60 = qJD(2) ^ 2;
t71 = t50 ^ 2 * t60;
t64 = t71 / 0.2e1;
t3 = t13 * t73 - t54 * t8;
t63 = qJD(2) * t70;
t11 = -t20 * t55 + t24 * t74;
t7 = pkin(4) * t44 - t11;
t61 = qJD(1) ^ 2;
t32 = qJD(5) + t33;
t31 = -t50 * t41 + t46;
t30 = t32 ^ 2 / 0.2e1;
t29 = t35 * t73 - t44 * t54;
t27 = t35 * t54 + t44 * t73;
t26 = t29 ^ 2 / 0.2e1;
t25 = t27 ^ 2 / 0.2e1;
t16 = t29 * t32;
t15 = t27 * t32;
t14 = t29 * t27;
t5 = pkin(5) * t27 + qJD(6) + t7;
t2 = -qJ(6) * t27 + t4;
t1 = pkin(5) * t32 - qJ(6) * t29 + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t61 / 0.2e1, 0, 0, 0, 0, 0, t60 / 0.2e1, t59 * t63, -t57 * t63, 0 (t53 ^ 2 / 0.2e1 + (t57 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1) * t51 ^ 2) * t61, t56 ^ 2 * t64, t56 * t58 * t71, t47 * t66, t58 ^ 2 * t64, t47 * t65, t47 ^ 2 / 0.2e1, t21 * t47 - t31 * t65, -t22 * t47 + t31 * t66 (-t21 * t56 + t22 * t58) * t68, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, -t35 * t44, t33 ^ 2 / 0.2e1, t33 * t44, t44 ^ 2 / 0.2e1, -t11 * t44 + t19 * t33, t12 * t44 + t19 * t35, -t11 * t35 - t12 * t33, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t26, -t14, t16, t25, -t15, t30, t27 * t7 + t3 * t32, t29 * t7 - t32 * t4, -t27 * t4 - t29 * t3, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t26, -t14, t16, t25, -t15, t30, t1 * t32 + t27 * t5, -t2 * t32 + t29 * t5, -t1 * t29 - t2 * t27, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
