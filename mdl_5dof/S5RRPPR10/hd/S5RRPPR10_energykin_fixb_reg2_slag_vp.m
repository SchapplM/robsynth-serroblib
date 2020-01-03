% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:39
% EndTime: 2019-12-31 19:44:39
% DurationCPUTime: 0.19s
% Computational Cost: add. (286->50), mult. (710->110), div. (0->0), fcn. (415->6), ass. (0->42)
t36 = qJD(1) ^ 2;
t52 = t36 / 0.2e1;
t51 = cos(qJ(5));
t32 = sin(pkin(8));
t34 = sin(qJ(2));
t47 = qJD(1) * t34;
t48 = cos(pkin(8));
t18 = -t48 * qJD(2) + t32 * t47;
t20 = t32 * qJD(2) + t48 * t47;
t50 = t20 * t18;
t35 = cos(qJ(2));
t49 = t35 * t36;
t17 = (-pkin(2) * t35 - qJ(3) * t34 - pkin(1)) * qJD(1);
t46 = t35 * qJD(1);
t24 = pkin(6) * t46 + qJD(2) * qJ(3);
t13 = t32 * t17 + t48 * t24;
t45 = t18 ^ 2 / 0.2e1;
t44 = qJD(1) * qJD(2);
t43 = t18 * t46;
t42 = t20 * t46;
t41 = t34 * t44;
t40 = t35 * t44;
t39 = qJD(2) * pkin(2) - pkin(6) * t47 - qJD(3);
t12 = t48 * t17 - t32 * t24;
t8 = -qJ(4) * t46 + t13;
t7 = pkin(3) * t46 + qJD(4) - t12;
t38 = t20 * qJ(4) + t39;
t33 = sin(qJ(5));
t31 = t35 ^ 2;
t30 = t34 ^ 2;
t26 = t31 * t52;
t25 = qJD(5) + t46;
t16 = t20 ^ 2 / 0.2e1;
t11 = t33 * t18 + t51 * t20;
t9 = -t51 * t18 + t33 * t20;
t6 = t18 * pkin(3) - t38;
t5 = t18 * pkin(7) + t8;
t4 = (-pkin(3) - pkin(4)) * t18 + t38;
t3 = pkin(4) * t46 - t20 * pkin(7) + t7;
t2 = t33 * t3 + t51 * t5;
t1 = t51 * t3 - t33 * t5;
t10 = [0, 0, 0, 0, 0, t52, 0, 0, 0, 0, t30 * t52, t34 * t49, t41, t26, t40, qJD(2) ^ 2 / 0.2e1, pkin(1) * t49 - pkin(6) * t41, -t36 * pkin(1) * t34 - pkin(6) * t40, (t30 + t31) * t36 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t31 / 0.2e1 + t30 / 0.2e1) * pkin(6) ^ 2) * t36, t16, -t50, -t42, t45, t43, t26, -t12 * t46 - t18 * t39, t13 * t46 - t20 * t39, -t12 * t20 - t13 * t18, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, t16, -t42, t50, t26, -t43, t45, t6 * t18 + t7 * t46, -t8 * t18 + t7 * t20, -t6 * t20 - t8 * t46, t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t9, t11 * t25, t9 ^ 2 / 0.2e1, -t9 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t4 * t9, t4 * t11 - t2 * t25, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg = t10;
