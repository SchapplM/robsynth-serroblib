% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR14_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:30
% EndTime: 2024-09-27 18:43:30
% DurationCPUTime: 0.09s
% Computational Cost: add. (648->46), mult. (1108->125), div. (0->0), fcn. (752->10), ass. (0->48)
t41 = cos(qJ(3));
t34 = qJD(1) + qJD(2);
t42 = cos(qJ(2));
t51 = pkin(1) * qJD(1);
t46 = t42 * t51;
t29 = t34 * pkin(2) + t46;
t36 = cos(pkin(5));
t55 = t29 * t36;
t24 = t41 * t55;
t40 = sin(qJ(2));
t47 = t40 * t51;
t35 = sin(pkin(5));
t53 = t34 * t35;
t26 = pkin(8) * t53 + t47;
t31 = t36 * t34 + qJD(3);
t39 = sin(qJ(3));
t12 = t31 * pkin(3) + t24 + (-pkin(9) * t53 - t26) * t39;
t17 = t41 * t26 + t39 * t55;
t52 = t34 * t41;
t48 = t35 * t52;
t15 = pkin(9) * t48 + t17;
t38 = sin(qJ(4));
t57 = cos(qJ(4));
t6 = t38 * t12 + t57 * t15;
t56 = cos(qJ(5));
t32 = t34 ^ 2;
t33 = t35 ^ 2;
t54 = t32 * t33;
t50 = t29 * t33 * t34;
t49 = t39 * t53;
t45 = t54 / 0.2e1;
t5 = t57 * t12 - t38 * t15;
t30 = qJD(4) + t31;
t19 = (-pkin(3) * t52 - t29) * t35;
t43 = qJD(1) ^ 2;
t37 = sin(qJ(5));
t28 = qJD(5) + t30;
t22 = (t38 * t41 + t57 * t39) * t53;
t20 = t38 * t49 - t57 * t48;
t16 = -t39 * t26 + t24;
t13 = t20 * pkin(4) + t19;
t11 = -t37 * t20 + t56 * t22;
t9 = t56 * t20 + t37 * t22;
t4 = -t20 * pkin(10) + t6;
t3 = t30 * pkin(4) - t22 * pkin(10) + t5;
t2 = t37 * t3 + t56 * t4;
t1 = t56 * t3 - t37 * t4;
t7 = [0, 0, 0, 0, 0, t43 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 / 0.2e1, t34 * t46, -t34 * t47, 0, (t40 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t43, t39 ^ 2 * t45, t39 * t41 * t54, t31 * t49, t41 ^ 2 * t45, t31 * t48, t31 ^ 2 / 0.2e1, t16 * t31 + t41 * t50, -t17 * t31 - t39 * t50, (-t16 * t39 + t17 * t41) * t53, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t33 * t29 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t30, t20 ^ 2 / 0.2e1, -t20 * t30, t30 ^ 2 / 0.2e1, t19 * t20 + t5 * t30, t19 * t22 - t6 * t30, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t9, t11 * t28, t9 ^ 2 / 0.2e1, -t9 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t13 * t9, t13 * t11 - t2 * t28, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg = t7;
