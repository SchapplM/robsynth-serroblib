% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:25
% EndTime: 2019-12-31 21:20:25
% DurationCPUTime: 0.19s
% Computational Cost: add. (309->49), mult. (743->110), div. (0->0), fcn. (456->6), ass. (0->46)
t32 = sin(qJ(3));
t33 = sin(qJ(2));
t34 = cos(qJ(3));
t35 = cos(qJ(2));
t19 = (t32 * t35 + t33 * t34) * qJD(1);
t36 = qJD(1) ^ 2;
t55 = t36 / 0.2e1;
t54 = pkin(3) + pkin(8);
t53 = -pkin(7) - pkin(6);
t52 = cos(qJ(5));
t46 = qJD(1) * t35;
t47 = qJD(1) * t33;
t17 = t32 * t47 - t34 * t46;
t51 = t19 * t17;
t28 = qJD(2) + qJD(3);
t50 = t19 * t28;
t49 = t28 * t17;
t48 = t35 * t36;
t23 = qJD(2) * pkin(2) + t53 * t47;
t24 = t53 * t46;
t10 = t32 * t23 - t34 * t24;
t45 = t17 ^ 2 / 0.2e1;
t44 = t19 ^ 2 / 0.2e1;
t43 = qJD(1) * qJD(2);
t42 = t33 * t43;
t41 = t35 * t43;
t9 = t34 * t23 + t32 * t24;
t8 = -t28 * qJ(4) - t10;
t40 = qJD(4) - t9;
t25 = (-pkin(2) * t35 - pkin(1)) * qJD(1);
t38 = -t19 * qJ(4) + t25;
t31 = sin(qJ(5));
t30 = t35 ^ 2;
t29 = t33 ^ 2;
t27 = t28 ^ 2 / 0.2e1;
t16 = qJD(5) + t19;
t13 = t31 * t17 + t52 * t28;
t11 = -t52 * t17 + t31 * t28;
t7 = -t28 * pkin(3) + t40;
t6 = t17 * pkin(3) + t38;
t5 = -t17 * pkin(4) - t8;
t4 = t54 * t17 + t38;
t3 = t19 * pkin(4) - t54 * t28 + t40;
t2 = t31 * t3 + t52 * t4;
t1 = t52 * t3 - t31 * t4;
t12 = [0, 0, 0, 0, 0, t55, 0, 0, 0, 0, t29 * t55, t33 * t48, t42, t30 * t55, t41, qJD(2) ^ 2 / 0.2e1, pkin(1) * t48 - pkin(6) * t42, -t36 * pkin(1) * t33 - pkin(6) * t41, (t29 + t30) * t36 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t30 / 0.2e1 + t29 / 0.2e1) * pkin(6) ^ 2) * t36, t44, -t51, t50, t45, -t49, t27, t25 * t17 + t9 * t28, -t10 * t28 + t25 * t19, -t10 * t17 - t9 * t19, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t27, -t50, t49, t44, -t51, t45, t8 * t17 + t7 * t19, -t6 * t17 + t7 * t28, -t6 * t19 - t8 * t28, t6 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t16, t11 ^ 2 / 0.2e1, -t11 * t16, t16 ^ 2 / 0.2e1, t1 * t16 + t5 * t11, t5 * t13 - t2 * t16, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t12;
