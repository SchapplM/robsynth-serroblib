% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:58
% EndTime: 2019-12-31 19:47:58
% DurationCPUTime: 0.18s
% Computational Cost: add. (286->50), mult. (668->113), div. (0->0), fcn. (357->6), ass. (0->42)
t40 = qJD(1) ^ 2;
t52 = t40 / 0.2e1;
t51 = cos(qJ(5));
t39 = cos(qJ(2));
t50 = t39 * t40;
t49 = -pkin(2) - qJ(4);
t38 = sin(qJ(2));
t42 = -qJ(3) * t38 - pkin(1);
t14 = (t49 * t39 + t42) * qJD(1);
t47 = t38 * qJD(1);
t46 = pkin(6) * t47 + qJD(3);
t15 = pkin(3) * t47 + t49 * qJD(2) + t46;
t35 = sin(pkin(8));
t36 = cos(pkin(8));
t6 = t36 * t14 + t35 * t15;
t48 = qJD(1) * t39;
t23 = -pkin(6) * t48 - qJD(2) * qJ(3);
t45 = qJD(1) * qJD(2);
t44 = t38 * t45;
t43 = t39 * t45;
t5 = -t35 * t14 + t36 * t15;
t17 = pkin(3) * t48 + qJD(4) - t23;
t37 = sin(qJ(5));
t34 = t39 ^ 2;
t33 = t38 ^ 2;
t31 = qJD(2) ^ 2 / 0.2e1;
t27 = t34 * t52;
t26 = t33 * t52;
t25 = qJD(5) + t47;
t24 = t38 * t50;
t22 = -qJD(2) * pkin(2) + t46;
t21 = t36 * qJD(2) - t35 * t48;
t19 = t35 * qJD(2) + t36 * t48;
t18 = (-pkin(2) * t39 + t42) * qJD(1);
t10 = t19 * pkin(4) + t17;
t9 = -t37 * t19 + t51 * t21;
t7 = t51 * t19 + t37 * t21;
t4 = -t19 * pkin(7) + t6;
t3 = pkin(4) * t47 - t21 * pkin(7) + t5;
t2 = t37 * t3 + t51 * t4;
t1 = t51 * t3 - t37 * t4;
t8 = [0, 0, 0, 0, 0, t52, 0, 0, 0, 0, t26, t24, t44, t27, t43, t31, pkin(1) * t50 - pkin(6) * t44, -t40 * pkin(1) * t38 - pkin(6) * t43, (t33 + t34) * t40 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t34 / 0.2e1 + t33 / 0.2e1) * pkin(6) ^ 2) * t40, t31, -t44, -t43, t26, t24, t27, (t22 * t38 - t23 * t39) * qJD(1), t22 * qJD(2) + t18 * t48, -t23 * qJD(2) - t18 * t47, t18 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * t47, t19 ^ 2 / 0.2e1, -t19 * t47, t26, t17 * t19 + t5 * t47, t17 * t21 - t6 * t47, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t25, t7 ^ 2 / 0.2e1, -t7 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t10 * t7, t10 * t9 - t2 * t25, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
