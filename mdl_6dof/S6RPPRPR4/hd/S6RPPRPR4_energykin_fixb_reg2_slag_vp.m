% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:07
% EndTime: 2019-03-09 01:47:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (421->52), mult. (804->120), div. (0->0), fcn. (435->8), ass. (0->41)
t48 = qJD(1) ^ 2;
t39 = t48 / 0.2e1;
t29 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t41 = sin(pkin(9));
t43 = cos(pkin(9));
t51 = qJ(2) * qJD(1);
t22 = t41 * t29 + t43 * t51;
t20 = -qJD(1) * pkin(7) + t22;
t47 = cos(qJ(4));
t37 = t47 * qJD(3);
t46 = sin(qJ(4));
t50 = qJ(5) * qJD(1);
t10 = qJD(4) * pkin(4) + t37 + (-t20 + t50) * t46;
t13 = t46 * qJD(3) + t47 * t20;
t11 = -t47 * t50 + t13;
t40 = sin(pkin(10));
t42 = cos(pkin(10));
t6 = t40 * t10 + t42 * t11;
t54 = cos(qJ(6));
t53 = qJD(1) * t46;
t52 = qJD(1) * t47;
t49 = qJD(1) * qJD(4);
t21 = t43 * t29 - t41 * t51;
t19 = qJD(1) * pkin(3) - t21;
t25 = t40 * t53 - t42 * t52;
t5 = t42 * t10 - t40 * t11;
t14 = pkin(4) * t52 + qJD(5) + t19;
t45 = sin(qJ(6));
t38 = qJD(4) ^ 2 / 0.2e1;
t34 = -qJD(1) * pkin(1) + qJD(2);
t26 = (-t40 * t47 - t42 * t46) * qJD(1);
t23 = -qJD(6) + t25;
t18 = t45 * qJD(4) + t54 * t26;
t16 = -t54 * qJD(4) + t45 * t26;
t12 = -t46 * t20 + t37;
t7 = -t25 * pkin(5) - t26 * pkin(8) + t14;
t4 = qJD(4) * pkin(8) + t6;
t3 = -qJD(4) * pkin(5) - t5;
t2 = t54 * t4 + t45 * t7;
t1 = -t45 * t4 + t54 * t7;
t8 = [0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, -t34 * qJD(1), 0, t48 * qJ(2), qJ(2) ^ 2 * t39 + t34 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t39, -t21 * qJD(1), t22 * qJD(1), 0, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t46 ^ 2 * t39, t46 * t48 * t47, -t46 * t49, t47 ^ 2 * t39, -t47 * t49, t38, t12 * qJD(4) + t19 * t52, -t13 * qJD(4) - t19 * t53 (t12 * t46 - t13 * t47) * qJD(1), t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, t26 * t25, t26 * qJD(4), t25 ^ 2 / 0.2e1, t25 * qJD(4), t38, t5 * qJD(4) - t14 * t25, -t6 * qJD(4) + t14 * t26, t6 * t25 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, -t18 * t23, t16 ^ 2 / 0.2e1, t16 * t23, t23 ^ 2 / 0.2e1, -t1 * t23 + t3 * t16, t3 * t18 + t2 * t23, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
