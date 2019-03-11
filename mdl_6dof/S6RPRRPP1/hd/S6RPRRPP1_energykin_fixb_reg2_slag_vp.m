% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:29:50
% EndTime: 2019-03-09 04:29:50
% DurationCPUTime: 0.16s
% Computational Cost: add. (501->56), mult. (1102->126), div. (0->0), fcn. (687->8), ass. (0->44)
t46 = qJD(1) ^ 2;
t39 = t46 / 0.2e1;
t40 = sin(pkin(10));
t53 = cos(pkin(10));
t41 = sin(pkin(9));
t32 = (pkin(1) * t41 + pkin(7)) * qJD(1);
t44 = sin(qJ(3));
t45 = cos(qJ(3));
t25 = t44 * qJD(2) + t45 * t32;
t22 = qJD(3) * pkin(8) + t25;
t42 = cos(pkin(9));
t48 = -pkin(1) * t42 - pkin(2);
t23 = (-pkin(3) * t45 - pkin(8) * t44 + t48) * qJD(1);
t43 = sin(qJ(4));
t56 = cos(qJ(4));
t10 = -t43 * t22 + t56 * t23;
t52 = qJD(1) * t44;
t31 = t43 * qJD(3) + t56 * t52;
t51 = t45 * qJD(1);
t35 = -qJD(4) + t51;
t7 = -t35 * pkin(4) - t31 * qJ(5) + t10;
t11 = t56 * t22 + t43 * t23;
t29 = -t56 * qJD(3) + t43 * t52;
t9 = -t29 * qJ(5) + t11;
t4 = t40 * t7 + t53 * t9;
t57 = pkin(1) * t46;
t15 = t53 * t29 + t40 * t31;
t17 = -t40 * t29 + t53 * t31;
t55 = t17 * t15;
t54 = t35 * t15;
t50 = t15 ^ 2 / 0.2e1;
t49 = qJD(1) * qJD(3);
t24 = t45 * qJD(2) - t44 * t32;
t3 = -t40 * t9 + t53 * t7;
t21 = -qJD(3) * pkin(3) - t24;
t13 = t29 * pkin(4) + qJD(5) + t21;
t34 = t35 ^ 2 / 0.2e1;
t33 = t48 * qJD(1);
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t35;
t5 = t15 * pkin(5) - t17 * qJ(6) + t13;
t2 = -t35 * qJ(6) + t4;
t1 = t35 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t42 * t57, -t41 * t57, 0, qJD(2) ^ 2 / 0.2e1 + (t41 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t46, t44 ^ 2 * t39, t45 * t46 * t44, t44 * t49, t45 ^ 2 * t39, t45 * t49, qJD(3) ^ 2 / 0.2e1, t24 * qJD(3) - t33 * t51, -t25 * qJD(3) + t33 * t52 (-t24 * t44 + t25 * t45) * qJD(1), t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t35, t29 ^ 2 / 0.2e1, t29 * t35, t34, -t10 * t35 + t21 * t29, t11 * t35 + t21 * t31, -t10 * t31 - t11 * t29, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14, -t55, -t12, t50, t54, t34, t13 * t15 - t3 * t35, t13 * t17 + t4 * t35, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t14, -t12, t55, t34, -t54, t50, t1 * t35 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 - t2 * t35, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
