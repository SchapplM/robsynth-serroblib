% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:21
% EndTime: 2019-03-09 00:45:21
% DurationCPUTime: 0.17s
% Computational Cost: add. (682->62), mult. (1509->154), div. (0->0), fcn. (1124->12), ass. (0->52)
t55 = qJD(2) ^ 2;
t66 = t55 / 0.2e1;
t51 = sin(qJ(2));
t45 = sin(pkin(6));
t63 = qJD(1) * t45;
t37 = qJD(2) * pkin(8) + t51 * t63;
t53 = cos(qJ(3));
t46 = cos(pkin(6));
t62 = qJD(1) * t46;
t40 = t53 * t62;
t50 = sin(qJ(3));
t22 = qJD(3) * pkin(3) + t40 + (-pkin(9) * qJD(2) - t37) * t50;
t29 = t53 * t37 + t50 * t62;
t60 = qJD(2) * t53;
t24 = pkin(9) * t60 + t29;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t12 = t49 * t22 + t52 * t24;
t44 = qJD(3) + qJD(4);
t10 = t44 * pkin(10) + t12;
t54 = cos(qJ(2));
t58 = t54 * t63;
t31 = -t58 + (-pkin(3) * t53 - pkin(2)) * qJD(2);
t61 = qJD(2) * t50;
t33 = t49 * t61 - t52 * t60;
t35 = (t49 * t53 + t50 * t52) * qJD(2);
t18 = t33 * pkin(4) - t35 * pkin(10) + t31;
t48 = sin(qJ(5));
t65 = cos(qJ(5));
t6 = t65 * t10 + t48 * t18;
t64 = cos(qJ(6));
t59 = qJD(2) * qJD(3);
t57 = qJD(2) * t63;
t5 = -t48 * t10 + t65 * t18;
t11 = t52 * t22 - t49 * t24;
t9 = -t44 * pkin(4) - t11;
t32 = qJD(5) + t33;
t56 = qJD(1) ^ 2;
t47 = sin(qJ(6));
t38 = -qJD(2) * pkin(2) - t58;
t30 = qJD(6) + t32;
t28 = -t50 * t37 + t40;
t27 = t65 * t35 + t48 * t44;
t25 = t48 * t35 - t65 * t44;
t15 = -t47 * t25 + t64 * t27;
t13 = t64 * t25 + t47 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(11) + t6;
t3 = t32 * pkin(5) - t27 * pkin(11) + t5;
t2 = t47 * t3 + t64 * t4;
t1 = t64 * t3 - t47 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t56 / 0.2e1, 0, 0, 0, 0, 0, t66, t54 * t57, -t51 * t57, 0 (t46 ^ 2 / 0.2e1 + (t51 ^ 2 / 0.2e1 + t54 ^ 2 / 0.2e1) * t45 ^ 2) * t56, t50 ^ 2 * t66, t50 * t55 * t53, t50 * t59, t53 ^ 2 * t66, t53 * t59, qJD(3) ^ 2 / 0.2e1, t28 * qJD(3) - t38 * t60, -t29 * qJD(3) + t38 * t61 (-t28 * t50 + t29 * t53) * qJD(2), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * t44, t33 ^ 2 / 0.2e1, -t33 * t44, t44 ^ 2 / 0.2e1, t11 * t44 + t31 * t33, -t12 * t44 + t31 * t35, -t11 * t35 - t12 * t33, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t32, t25 ^ 2 / 0.2e1, -t25 * t32, t32 ^ 2 / 0.2e1, t9 * t25 + t5 * t32, t9 * t27 - t6 * t32, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t30, t13 ^ 2 / 0.2e1, -t13 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t13, t7 * t15 - t2 * t30, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
