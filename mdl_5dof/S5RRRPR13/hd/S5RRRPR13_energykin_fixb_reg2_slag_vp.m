% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR13_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:31
% EndTime: 2019-12-31 21:46:32
% DurationCPUTime: 0.17s
% Computational Cost: add. (454->53), mult. (1141->117), div. (0->0), fcn. (815->8), ass. (0->48)
t61 = pkin(3) + pkin(9);
t60 = cos(qJ(5));
t54 = cos(pkin(5)) * qJD(1);
t35 = qJD(2) + t54;
t40 = sin(qJ(3));
t42 = cos(qJ(3));
t41 = sin(qJ(2));
t37 = sin(pkin(5));
t55 = qJD(1) * t37;
t49 = t41 * t55;
t22 = -t42 * t35 + t40 * t49;
t24 = t40 * t35 + t42 * t49;
t59 = t24 * t22;
t43 = cos(qJ(2));
t48 = t43 * t55;
t30 = -qJD(3) + t48;
t58 = t24 * t30;
t57 = t30 * t22;
t44 = qJD(1) ^ 2;
t56 = t37 ^ 2 * t44;
t50 = pkin(1) * t54;
t27 = pkin(7) * t48 + t41 * t50;
t18 = t35 * pkin(8) + t27;
t20 = (-pkin(2) * t43 - pkin(8) * t41 - pkin(1)) * t55;
t10 = t42 * t18 + t40 * t20;
t53 = t22 ^ 2 / 0.2e1;
t52 = t24 ^ 2 / 0.2e1;
t51 = t43 * t56;
t47 = t56 / 0.2e1;
t9 = -t40 * t18 + t42 * t20;
t8 = t30 * qJ(4) - t10;
t46 = qJD(4) - t9;
t26 = -pkin(7) * t49 + t43 * t50;
t17 = -t35 * pkin(2) - t26;
t45 = -t24 * qJ(4) + t17;
t39 = sin(qJ(5));
t28 = t30 ^ 2 / 0.2e1;
t21 = qJD(5) + t24;
t13 = t39 * t22 - t60 * t30;
t11 = -t60 * t22 - t39 * t30;
t7 = t30 * pkin(3) + t46;
t6 = t22 * pkin(3) + t45;
t5 = -t22 * pkin(4) - t8;
t4 = t61 * t22 + t45;
t3 = t24 * pkin(4) + t61 * t30 + t46;
t2 = t39 * t3 + t60 * t4;
t1 = t60 * t3 - t39 * t4;
t12 = [0, 0, 0, 0, 0, t44 / 0.2e1, 0, 0, 0, 0, t41 ^ 2 * t47, t41 * t51, t35 * t49, t43 ^ 2 * t47, t35 * t48, t35 ^ 2 / 0.2e1, pkin(1) * t51 + t26 * t35, -pkin(1) * t41 * t56 - t27 * t35, (-t26 * t41 + t27 * t43) * t55, t27 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t47, t52, -t59, -t58, t53, t57, t28, t17 * t22 - t9 * t30, t10 * t30 + t17 * t24, -t10 * t22 - t9 * t24, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t28, t58, -t57, t52, -t59, t53, t8 * t22 + t7 * t24, -t6 * t22 - t7 * t30, -t6 * t24 + t8 * t30, t6 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t21, t11 ^ 2 / 0.2e1, -t11 * t21, t21 ^ 2 / 0.2e1, t1 * t21 + t5 * t11, t5 * t13 - t2 * t21, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t12;
