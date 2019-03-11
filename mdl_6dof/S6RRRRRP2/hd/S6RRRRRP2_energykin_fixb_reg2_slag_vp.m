% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:03:09
% EndTime: 2019-03-10 01:03:09
% DurationCPUTime: 0.22s
% Computational Cost: add. (909->62), mult. (2104->143), div. (0->0), fcn. (1532->8), ass. (0->53)
t52 = qJD(1) ^ 2;
t67 = t52 / 0.2e1;
t66 = -pkin(8) - pkin(7);
t49 = sin(qJ(3));
t51 = cos(qJ(2));
t58 = qJD(1) * t51;
t50 = sin(qJ(2));
t59 = qJD(1) * t50;
t65 = cos(qJ(3));
t34 = t49 * t59 - t65 * t58;
t36 = (t49 * t51 + t65 * t50) * qJD(1);
t48 = sin(qJ(4));
t64 = cos(qJ(4));
t25 = t64 * t34 + t48 * t36;
t27 = -t48 * t34 + t64 * t36;
t40 = (-pkin(2) * t51 - pkin(1)) * qJD(1);
t30 = t34 * pkin(3) + t40;
t12 = t25 * pkin(4) - t27 * pkin(10) + t30;
t47 = sin(qJ(5));
t63 = cos(qJ(5));
t38 = qJD(2) * pkin(2) + t66 * t59;
t39 = t66 * t58;
t28 = t65 * t38 + t49 * t39;
t44 = qJD(2) + qJD(3);
t15 = t44 * pkin(3) - t36 * pkin(9) + t28;
t29 = t49 * t38 - t65 * t39;
t19 = -t34 * pkin(9) + t29;
t10 = t48 * t15 + t64 * t19;
t43 = qJD(4) + t44;
t8 = t43 * pkin(10) + t10;
t4 = t47 * t12 + t63 * t8;
t20 = t47 * t27 - t63 * t43;
t22 = t63 * t27 + t47 * t43;
t62 = t22 * t20;
t24 = qJD(5) + t25;
t61 = t24 * t20;
t60 = t51 * t52;
t57 = t20 ^ 2 / 0.2e1;
t56 = qJD(1) * qJD(2);
t55 = t50 * t56;
t54 = t51 * t56;
t9 = t64 * t15 - t48 * t19;
t3 = t63 * t12 - t47 * t8;
t7 = -t43 * pkin(4) - t9;
t46 = t51 ^ 2;
t45 = t50 ^ 2;
t23 = t24 ^ 2 / 0.2e1;
t18 = t22 ^ 2 / 0.2e1;
t13 = t22 * t24;
t5 = t20 * pkin(5) - t22 * qJ(6) + t7;
t2 = t24 * qJ(6) + t4;
t1 = -t24 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t45 * t67, t50 * t60, t55, t46 * t67, t54, qJD(2) ^ 2 / 0.2e1, pkin(1) * t60 - pkin(7) * t55, -t52 * pkin(1) * t50 - pkin(7) * t54 (t45 + t46) * t52 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t46 / 0.2e1 + t45 / 0.2e1) * pkin(7) ^ 2) * t52, t36 ^ 2 / 0.2e1, -t36 * t34, t36 * t44, t34 ^ 2 / 0.2e1, -t34 * t44, t44 ^ 2 / 0.2e1, t28 * t44 + t40 * t34, -t29 * t44 + t40 * t36, -t28 * t36 - t29 * t34, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t43, t25 ^ 2 / 0.2e1, -t25 * t43, t43 ^ 2 / 0.2e1, t30 * t25 + t9 * t43, -t10 * t43 + t30 * t27, -t10 * t25 - t9 * t27, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t18, -t62, t13, t57, -t61, t23, t7 * t20 + t3 * t24, t7 * t22 - t4 * t24, -t4 * t20 - t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t18, t13, t62, t23, t61, t57, -t1 * t24 + t5 * t20, t1 * t22 - t2 * t20, t2 * t24 - t5 * t22, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
