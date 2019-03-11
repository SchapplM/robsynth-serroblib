% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:23:43
% EndTime: 2019-03-09 08:23:43
% DurationCPUTime: 0.18s
% Computational Cost: add. (418->64), mult. (986->126), div. (0->0), fcn. (588->6), ass. (0->49)
t39 = sin(pkin(9));
t41 = sin(qJ(2));
t54 = qJD(1) * t41;
t55 = cos(pkin(9));
t22 = -t55 * qJD(2) + t39 * t54;
t63 = (pkin(3) + qJ(5)) * t22;
t43 = qJD(1) ^ 2;
t62 = t43 / 0.2e1;
t61 = pkin(4) + pkin(8);
t60 = cos(qJ(6));
t24 = t39 * qJD(2) + t55 * t54;
t59 = t22 * t24;
t42 = cos(qJ(2));
t58 = t42 * t43;
t56 = -pkin(5) - qJ(4);
t21 = (-pkin(2) * t42 - qJ(3) * t41 - pkin(1)) * qJD(1);
t53 = t42 * qJD(1);
t29 = pkin(7) * t53 + qJD(2) * qJ(3);
t16 = t39 * t21 + t55 * t29;
t19 = t22 ^ 2 / 0.2e1;
t20 = t24 ^ 2 / 0.2e1;
t52 = qJD(1) * qJD(2);
t51 = t22 * t53;
t50 = t24 * t53;
t49 = t41 * t52;
t48 = t42 * t52;
t47 = qJD(2) * pkin(2) - pkin(7) * t54 - qJD(3);
t15 = t55 * t21 - t39 * t29;
t11 = qJ(4) * t53 - t16;
t10 = pkin(3) * t53 + qJD(4) - t15;
t46 = t24 * qJ(4) + t47;
t45 = qJ(5) * t53 + t10;
t40 = sin(qJ(6));
t38 = t42 ^ 2;
t37 = t41 ^ 2;
t32 = t38 * t62;
t30 = -qJD(6) + t53;
t14 = t60 * t22 + t40 * t24;
t12 = t40 * t22 - t60 * t24;
t9 = t22 * pkin(3) - t46;
t8 = -t22 * pkin(4) + qJD(5) - t11;
t7 = t46 - t63;
t6 = t24 * pkin(4) + t45;
t5 = t56 * t24 - t47 + t63;
t4 = -t61 * t22 + t56 * t53 + qJD(5) + t16;
t3 = t61 * t24 + t45;
t2 = t60 * t3 + t40 * t4;
t1 = -t40 * t3 + t60 * t4;
t13 = [0, 0, 0, 0, 0, t62, 0, 0, 0, 0, t37 * t62, t41 * t58, t49, t32, t48, qJD(2) ^ 2 / 0.2e1, pkin(1) * t58 - pkin(7) * t49, -t43 * pkin(1) * t41 - pkin(7) * t48 (t37 + t38) * t43 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t38 / 0.2e1 + t37 / 0.2e1) * pkin(7) ^ 2) * t43, t20, -t59, -t50, t19, t51, t32, -t15 * t53 - t22 * t47, t16 * t53 - t24 * t47, -t15 * t24 - t16 * t22, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1, t32, t50, -t51, t20, -t59, t19, t10 * t24 + t11 * t22, -t10 * t53 - t9 * t22, t11 * t53 - t9 * t24, t9 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t19, t51, t59, t32, t50, t20, t7 * t24 - t8 * t53, t8 * t22 - t6 * t24, -t7 * t22 + t6 * t53, t6 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, -t14 * t30, t12 ^ 2 / 0.2e1, t12 * t30, t30 ^ 2 / 0.2e1, -t1 * t30 + t5 * t12, t5 * t14 + t2 * t30, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t13;
