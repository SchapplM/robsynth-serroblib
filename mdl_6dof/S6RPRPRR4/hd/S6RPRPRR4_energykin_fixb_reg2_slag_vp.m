% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:26
% EndTime: 2019-03-09 03:46:26
% DurationCPUTime: 0.15s
% Computational Cost: add. (395->57), mult. (840->129), div. (0->0), fcn. (459->8), ass. (0->48)
t47 = qJD(1) ^ 2;
t38 = t47 / 0.2e1;
t59 = -pkin(3) - pkin(8);
t58 = pkin(1) * t47;
t40 = sin(pkin(10));
t28 = (pkin(1) * t40 + pkin(7)) * qJD(1);
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t20 = t46 * qJD(2) - t44 * t28;
t50 = qJD(4) - t20;
t55 = t44 * qJD(1);
t13 = pkin(4) * t55 + t59 * qJD(3) + t50;
t41 = cos(pkin(10));
t53 = -pkin(1) * t41 - pkin(2);
t49 = -qJ(4) * t44 + t53;
t16 = (t59 * t46 + t49) * qJD(1);
t43 = sin(qJ(5));
t45 = cos(qJ(5));
t6 = t43 * t13 + t45 * t16;
t57 = cos(qJ(6));
t21 = t44 * qJD(2) + t46 * t28;
t56 = qJD(1) * t46;
t54 = qJD(1) * qJD(3);
t18 = -qJD(3) * qJ(4) - t21;
t5 = t45 * t13 - t43 * t16;
t52 = t44 * t54;
t51 = t46 * t54;
t15 = pkin(4) * t56 - t18;
t32 = qJD(5) + t55;
t42 = sin(qJ(6));
t37 = qJD(3) ^ 2 / 0.2e1;
t34 = t46 ^ 2 * t38;
t33 = t44 ^ 2 * t38;
t31 = t46 * t47 * t44;
t30 = qJD(6) + t32;
t29 = t53 * qJD(1);
t27 = t45 * qJD(3) - t43 * t56;
t25 = t43 * qJD(3) + t45 * t56;
t19 = (-pkin(3) * t46 + t49) * qJD(1);
t17 = -qJD(3) * pkin(3) + t50;
t12 = -t42 * t25 + t57 * t27;
t10 = t57 * t25 + t42 * t27;
t7 = t25 * pkin(5) + t15;
t4 = -t25 * pkin(9) + t6;
t3 = t32 * pkin(5) - t27 * pkin(9) + t5;
t2 = t42 * t3 + t57 * t4;
t1 = t57 * t3 - t42 * t4;
t8 = [0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41 * t58, -t40 * t58, 0, qJD(2) ^ 2 / 0.2e1 + (t40 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t47, t33, t31, t52, t34, t51, t37, t20 * qJD(3) - t29 * t56, -t21 * qJD(3) + t29 * t55 (-t20 * t44 + t21 * t46) * qJD(1), t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t37, -t52, -t51, t33, t31, t34 (t17 * t44 - t18 * t46) * qJD(1), t17 * qJD(3) + t19 * t56, -t18 * qJD(3) - t19 * t55, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t32, t25 ^ 2 / 0.2e1, -t25 * t32, t32 ^ 2 / 0.2e1, t15 * t25 + t5 * t32, t15 * t27 - t6 * t32, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t12 ^ 2 / 0.2e1, -t12 * t10, t12 * t30, t10 ^ 2 / 0.2e1, -t10 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t10, t7 * t12 - t2 * t30, -t1 * t12 - t2 * t10, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
