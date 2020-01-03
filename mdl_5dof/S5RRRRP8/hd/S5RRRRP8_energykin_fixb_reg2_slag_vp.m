% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:07
% EndTime: 2019-12-31 22:02:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (388->49), mult. (905->111), div. (0->0), fcn. (579->6), ass. (0->42)
t46 = qJD(1) ^ 2;
t55 = t46 / 0.2e1;
t43 = sin(qJ(2));
t45 = cos(qJ(2));
t25 = (-pkin(2) * t45 - pkin(7) * t43 - pkin(1)) * qJD(1);
t51 = t45 * qJD(1);
t32 = pkin(6) * t51 + qJD(2) * pkin(7);
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t20 = t42 * t25 + t44 * t32;
t52 = qJD(1) * t43;
t26 = -t44 * qJD(2) + t42 * t52;
t11 = -t26 * pkin(8) + t20;
t41 = sin(qJ(4));
t54 = cos(qJ(4));
t19 = t44 * t25 - t42 * t32;
t28 = t42 * qJD(2) + t44 * t52;
t35 = -qJD(3) + t51;
t9 = -t35 * pkin(3) - t28 * pkin(8) + t19;
t4 = t54 * t11 + t41 * t9;
t53 = t45 * t46;
t50 = qJD(1) * qJD(2);
t3 = -t41 * t11 + t54 * t9;
t49 = t43 * t50;
t48 = t45 * t50;
t31 = -qJD(2) * pkin(2) + pkin(6) * t52;
t21 = t26 * pkin(3) + t31;
t40 = t45 ^ 2;
t39 = t43 ^ 2;
t33 = -qJD(4) + t35;
t30 = t33 ^ 2 / 0.2e1;
t18 = -t41 * t26 + t54 * t28;
t16 = t54 * t26 + t41 * t28;
t15 = t18 ^ 2 / 0.2e1;
t14 = t16 ^ 2 / 0.2e1;
t13 = t18 * t33;
t12 = t16 * t33;
t6 = t16 * pkin(4) + qJD(5) + t21;
t5 = t18 * t16;
t2 = -t16 * qJ(5) + t4;
t1 = -t33 * pkin(4) - t18 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, t55, 0, 0, 0, 0, t39 * t55, t43 * t53, t49, t40 * t55, t48, qJD(2) ^ 2 / 0.2e1, pkin(1) * t53 - pkin(6) * t49, -t46 * pkin(1) * t43 - pkin(6) * t48, (t39 + t40) * t46 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t40 / 0.2e1 + t39 / 0.2e1) * pkin(6) ^ 2) * t46, t28 ^ 2 / 0.2e1, -t28 * t26, -t28 * t35, t26 ^ 2 / 0.2e1, t26 * t35, t35 ^ 2 / 0.2e1, -t19 * t35 + t31 * t26, t20 * t35 + t31 * t28, -t19 * t28 - t20 * t26, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t15, -t5, -t13, t14, t12, t30, t21 * t16 - t3 * t33, t21 * t18 + t4 * t33, -t4 * t16 - t3 * t18, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t15, -t5, -t13, t14, t12, t30, -t1 * t33 + t6 * t16, t6 * t18 + t2 * t33, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg = t7;
