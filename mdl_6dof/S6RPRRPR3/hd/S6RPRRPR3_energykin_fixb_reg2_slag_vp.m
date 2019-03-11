% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:14
% EndTime: 2019-03-09 05:07:14
% DurationCPUTime: 0.17s
% Computational Cost: add. (420->57), mult. (896->128), div. (0->0), fcn. (526->8), ass. (0->48)
t47 = qJD(1) ^ 2;
t40 = t47 / 0.2e1;
t62 = pkin(4) + pkin(5);
t61 = pkin(1) * t47;
t60 = cos(qJ(4));
t59 = cos(qJ(6));
t44 = sin(qJ(4));
t45 = sin(qJ(3));
t56 = qJD(1) * t45;
t27 = -t60 * qJD(3) + t44 * t56;
t29 = t44 * qJD(3) + t60 * t56;
t58 = t29 * t27;
t46 = cos(qJ(3));
t55 = t46 * qJD(1);
t35 = -qJD(4) + t55;
t57 = t35 * t27;
t41 = sin(pkin(10));
t30 = (pkin(1) * t41 + pkin(7)) * qJD(1);
t22 = t45 * qJD(2) + t46 * t30;
t18 = qJD(3) * pkin(8) + t22;
t42 = cos(pkin(10));
t52 = -pkin(1) * t42 - pkin(2);
t19 = (-pkin(3) * t46 - pkin(8) * t45 + t52) * qJD(1);
t10 = t60 * t18 + t44 * t19;
t21 = t46 * qJD(2) - t45 * t30;
t54 = t27 ^ 2 / 0.2e1;
t53 = qJD(1) * qJD(3);
t7 = -t35 * qJ(5) + t10;
t51 = qJD(3) * pkin(3) + t21;
t9 = -t44 * t18 + t60 * t19;
t50 = qJD(5) - t9;
t49 = t29 * qJ(5) + t51;
t43 = sin(qJ(6));
t34 = qJD(6) + t35;
t32 = t35 ^ 2 / 0.2e1;
t31 = t52 * qJD(1);
t24 = t29 ^ 2 / 0.2e1;
t20 = t29 * t35;
t13 = t43 * t27 + t59 * t29;
t11 = -t59 * t27 + t43 * t29;
t8 = t27 * pkin(4) - t49;
t6 = t35 * pkin(4) + t50;
t5 = -t62 * t27 + t49;
t4 = t27 * pkin(9) + t7;
t3 = -t29 * pkin(9) + t62 * t35 + t50;
t2 = t43 * t3 + t59 * t4;
t1 = t59 * t3 - t43 * t4;
t12 = [0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t42 * t61, -t41 * t61, 0, qJD(2) ^ 2 / 0.2e1 + (t41 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t47, t45 ^ 2 * t40, t45 * t47 * t46, t45 * t53, t46 ^ 2 * t40, t46 * t53, qJD(3) ^ 2 / 0.2e1, t21 * qJD(3) - t31 * t55, -t22 * qJD(3) + t31 * t56 (-t21 * t45 + t22 * t46) * qJD(1), t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t24, -t58, -t20, t54, t57, t32, -t27 * t51 - t9 * t35, t10 * t35 - t29 * t51, -t10 * t27 - t9 * t29, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1, t24, -t20, t58, t32, -t57, t54, t8 * t27 + t6 * t35, -t7 * t27 + t6 * t29, -t8 * t29 - t7 * t35, t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t34, t11 ^ 2 / 0.2e1, -t11 * t34, t34 ^ 2 / 0.2e1, t1 * t34 + t5 * t11, t5 * t13 - t2 * t34, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
