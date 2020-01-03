% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:26
% EndTime: 2019-12-31 20:04:26
% DurationCPUTime: 0.14s
% Computational Cost: add. (196->46), mult. (489->95), div. (0->0), fcn. (242->4), ass. (0->40)
t44 = qJD(1) ^ 2;
t53 = t44 / 0.2e1;
t41 = sin(qJ(2));
t51 = qJD(1) * t41;
t49 = pkin(6) * t51 + qJD(3);
t12 = -pkin(7) * t51 + (-pkin(2) - pkin(3)) * qJD(2) + t49;
t43 = cos(qJ(2));
t50 = qJD(1) * t43;
t22 = pkin(6) * t50 + qJD(2) * qJ(3);
t19 = -pkin(7) * t50 + t22;
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t4 = t40 * t12 + t42 * t19;
t52 = t43 * t44;
t48 = qJD(1) * qJD(2);
t47 = t41 * t52;
t20 = -qJD(1) * pkin(1) - pkin(2) * t50 - qJ(3) * t51;
t3 = t42 * t12 - t40 * t19;
t26 = t41 * t48;
t46 = t43 * t48;
t9 = pkin(3) * t50 - t20;
t38 = t43 ^ 2;
t37 = t41 ^ 2;
t35 = qJD(2) ^ 2 / 0.2e1;
t33 = qJD(2) - qJD(4);
t28 = t33 ^ 2 / 0.2e1;
t25 = t38 * t53;
t24 = t37 * t53;
t21 = -qJD(2) * pkin(2) + t49;
t18 = -t40 * t50 + t42 * t51;
t16 = (-t40 * t41 - t42 * t43) * qJD(1);
t15 = t18 ^ 2 / 0.2e1;
t14 = t16 ^ 2 / 0.2e1;
t11 = t18 * t33;
t10 = t16 * t33;
t6 = t18 * t16;
t5 = -t16 * pkin(4) + qJD(5) + t9;
t2 = t16 * qJ(5) + t4;
t1 = -t33 * pkin(4) - t18 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, t53, 0, 0, 0, 0, t24, t47, t26, t25, t46, t35, pkin(1) * t52 - pkin(6) * t26, -t44 * pkin(1) * t41 - pkin(6) * t46, (t37 + t38) * t44 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t38 / 0.2e1 + t37 / 0.2e1) * pkin(6) ^ 2) * t44, t24, t26, -t47, t35, -t46, t25, -t21 * qJD(2) - t20 * t50, (t21 * t41 + t22 * t43) * qJD(1), t22 * qJD(2) - t20 * t51, t22 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t15, t6, -t11, t14, -t10, t28, -t9 * t16 - t3 * t33, t9 * t18 + t4 * t33, t4 * t16 - t3 * t18, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15, t6, -t11, t14, -t10, t28, -t1 * t33 - t5 * t16, t5 * t18 + t2 * t33, -t1 * t18 + t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
