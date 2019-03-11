% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:49:45
% EndTime: 2019-03-09 21:49:46
% DurationCPUTime: 0.20s
% Computational Cost: add. (841->63), mult. (1980->133), div. (0->0), fcn. (1489->8), ass. (0->51)
t52 = cos(qJ(2));
t51 = sin(qJ(2));
t47 = sin(pkin(6));
t62 = qJD(1) * t47;
t58 = t51 * t62;
t61 = cos(pkin(6)) * qJD(1);
t59 = pkin(1) * t61;
t36 = -pkin(8) * t58 + t52 * t59;
t45 = qJD(2) + t61;
t28 = -t45 * pkin(2) - t36;
t50 = sin(qJ(3));
t67 = cos(qJ(3));
t33 = -t67 * t45 + t50 * t58;
t35 = t50 * t45 + t67 * t58;
t10 = t33 * pkin(3) - t35 * pkin(10) + t28;
t57 = t52 * t62;
t37 = pkin(8) * t57 + t51 * t59;
t29 = t45 * pkin(9) + t37;
t31 = (-pkin(2) * t52 - pkin(9) * t51 - pkin(1)) * t62;
t18 = t67 * t29 + t50 * t31;
t40 = -qJD(3) + t57;
t14 = -t40 * pkin(10) + t18;
t49 = sin(qJ(4));
t66 = cos(qJ(4));
t8 = t49 * t10 + t66 * t14;
t21 = t49 * t35 + t66 * t40;
t23 = t66 * t35 - t49 * t40;
t65 = t21 * t23;
t32 = qJD(4) + t33;
t16 = t23 * t32;
t15 = t32 * t21;
t53 = qJD(1) ^ 2;
t64 = t47 ^ 2 * t53;
t63 = pkin(4) + qJ(6);
t19 = t21 ^ 2 / 0.2e1;
t20 = t23 ^ 2 / 0.2e1;
t60 = t52 * t64;
t56 = t64 / 0.2e1;
t5 = -t32 * qJ(5) - t8;
t7 = t66 * t10 - t49 * t14;
t17 = -t50 * t29 + t67 * t31;
t55 = qJD(5) - t7;
t13 = t40 * pkin(3) - t17;
t54 = -t23 * qJ(5) + t13;
t30 = t32 ^ 2 / 0.2e1;
t6 = t21 * pkin(4) + t54;
t4 = -t32 * pkin(4) + t55;
t3 = t63 * t21 + t54;
t2 = -t21 * pkin(5) + qJD(6) - t5;
t1 = t23 * pkin(5) - t63 * t32 + t55;
t9 = [0, 0, 0, 0, 0, t53 / 0.2e1, 0, 0, 0, 0, t51 ^ 2 * t56, t51 * t60, t45 * t58, t52 ^ 2 * t56, t45 * t57, t45 ^ 2 / 0.2e1, pkin(1) * t60 + t36 * t45, -pkin(1) * t51 * t64 - t37 * t45 (-t36 * t51 + t37 * t52) * t62, t37 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t56, t35 ^ 2 / 0.2e1, -t35 * t33, -t35 * t40, t33 ^ 2 / 0.2e1, t33 * t40, t40 ^ 2 / 0.2e1, -t17 * t40 + t28 * t33, t18 * t40 + t28 * t35, -t17 * t35 - t18 * t33, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t20, -t65, t16, t19, -t15, t30, t13 * t21 + t7 * t32, t13 * t23 - t8 * t32, -t8 * t21 - t7 * t23, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t30, -t16, t15, t20, -t65, t19, t5 * t21 + t4 * t23, -t6 * t21 + t4 * t32, -t6 * t23 - t5 * t32, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t30, t15, t16, t19, t65, t20, t1 * t23 - t2 * t21, t2 * t32 - t3 * t23, -t1 * t32 + t3 * t21, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg  = t9;
