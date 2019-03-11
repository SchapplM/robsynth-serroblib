% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:46:52
% EndTime: 2019-03-08 22:46:52
% DurationCPUTime: 0.15s
% Computational Cost: add. (501->57), mult. (1149->130), div. (0->0), fcn. (815->10), ass. (0->49)
t48 = qJD(2) ^ 2;
t62 = t48 / 0.2e1;
t40 = sin(pkin(11));
t58 = cos(pkin(11));
t45 = sin(qJ(2));
t41 = sin(pkin(6));
t57 = qJD(1) * t41;
t32 = qJD(2) * pkin(8) + t45 * t57;
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t42 = cos(pkin(6));
t56 = qJD(1) * t42;
t24 = t46 * t32 + t44 * t56;
t20 = qJD(3) * pkin(9) + t24;
t47 = cos(qJ(2));
t51 = t47 * t57;
t25 = -t51 + (-pkin(3) * t46 - pkin(9) * t44 - pkin(2)) * qJD(2);
t43 = sin(qJ(4));
t61 = cos(qJ(4));
t10 = -t43 * t20 + t61 * t25;
t55 = qJD(2) * t44;
t31 = t43 * qJD(3) + t61 * t55;
t54 = t46 * qJD(2);
t36 = -qJD(4) + t54;
t7 = -t36 * pkin(4) - t31 * qJ(5) + t10;
t11 = t61 * t20 + t43 * t25;
t29 = -t61 * qJD(3) + t43 * t55;
t9 = -t29 * qJ(5) + t11;
t4 = t40 * t7 + t58 * t9;
t15 = t58 * t29 + t40 * t31;
t17 = -t40 * t29 + t58 * t31;
t60 = t17 * t15;
t59 = t36 * t15;
t53 = t15 ^ 2 / 0.2e1;
t52 = qJD(2) * qJD(3);
t50 = qJD(2) * t57;
t23 = -t44 * t32 + t46 * t56;
t3 = -t40 * t9 + t58 * t7;
t19 = -qJD(3) * pkin(3) - t23;
t12 = t29 * pkin(4) + qJD(5) + t19;
t49 = qJD(1) ^ 2;
t34 = t36 ^ 2 / 0.2e1;
t33 = -qJD(2) * pkin(2) - t51;
t14 = t17 ^ 2 / 0.2e1;
t13 = t17 * t36;
t5 = t15 * pkin(5) - t17 * qJ(6) + t12;
t2 = -t36 * qJ(6) + t4;
t1 = t36 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t49 / 0.2e1, 0, 0, 0, 0, 0, t62, t47 * t50, -t45 * t50, 0 (t42 ^ 2 / 0.2e1 + (t45 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1) * t41 ^ 2) * t49, t44 ^ 2 * t62, t44 * t48 * t46, t44 * t52, t46 ^ 2 * t62, t46 * t52, qJD(3) ^ 2 / 0.2e1, t23 * qJD(3) - t33 * t54, -t24 * qJD(3) + t33 * t55 (-t23 * t44 + t24 * t46) * qJD(2), t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t36, t29 ^ 2 / 0.2e1, t29 * t36, t34, -t10 * t36 + t19 * t29, t11 * t36 + t19 * t31, -t10 * t31 - t11 * t29, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t14, -t60, -t13, t53, t59, t34, t12 * t15 - t3 * t36, t12 * t17 + t4 * t36, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t14, -t13, t60, t34, -t59, t53, t1 * t36 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 - t2 * t36, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
