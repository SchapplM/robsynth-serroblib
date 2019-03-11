% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:59:16
% EndTime: 2019-03-09 17:59:16
% DurationCPUTime: 0.20s
% Computational Cost: add. (766->63), mult. (1802->130), div. (0->0), fcn. (1329->8), ass. (0->57)
t73 = pkin(3) + pkin(10);
t48 = sin(qJ(5));
t64 = cos(pkin(6)) * qJD(1);
t44 = qJD(2) + t64;
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t50 = sin(qJ(2));
t46 = sin(pkin(6));
t65 = qJD(1) * t46;
t58 = t50 * t65;
t32 = t49 * t44 + t51 * t58;
t52 = cos(qJ(2));
t57 = t52 * t65;
t39 = -qJD(3) + t57;
t59 = pkin(1) * t64;
t35 = pkin(8) * t57 + t50 * t59;
t25 = t44 * pkin(9) + t35;
t28 = (-pkin(2) * t52 - pkin(9) * t50 - pkin(1)) * t65;
t15 = -t49 * t25 + t51 * t28;
t55 = qJD(4) - t15;
t7 = t32 * pkin(4) + t73 * t39 + t55;
t72 = cos(qJ(5));
t30 = -t51 * t44 + t49 * t58;
t34 = -pkin(8) * t58 + t52 * t59;
t24 = -t44 * pkin(2) - t34;
t54 = -t32 * qJ(4) + t24;
t9 = t73 * t30 + t54;
t4 = t48 * t7 + t72 * t9;
t18 = -t72 * t30 - t48 * t39;
t20 = t48 * t30 - t72 * t39;
t71 = t20 * t18;
t29 = qJD(5) + t32;
t70 = t29 * t18;
t69 = t32 * t30;
t68 = t32 * t39;
t67 = t39 * t30;
t53 = qJD(1) ^ 2;
t66 = t46 ^ 2 * t53;
t16 = t51 * t25 + t49 * t28;
t63 = t18 ^ 2 / 0.2e1;
t62 = t30 ^ 2 / 0.2e1;
t61 = t32 ^ 2 / 0.2e1;
t60 = t52 * t66;
t13 = t39 * qJ(4) - t16;
t56 = t66 / 0.2e1;
t10 = -t30 * pkin(4) - t13;
t3 = -t48 * t9 + t72 * t7;
t36 = t39 ^ 2 / 0.2e1;
t26 = t29 ^ 2 / 0.2e1;
t17 = t20 ^ 2 / 0.2e1;
t14 = t20 * t29;
t12 = t39 * pkin(3) + t55;
t11 = t30 * pkin(3) + t54;
t5 = t18 * pkin(5) - t20 * qJ(6) + t10;
t2 = t29 * qJ(6) + t4;
t1 = -t29 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t53 / 0.2e1, 0, 0, 0, 0, t50 ^ 2 * t56, t50 * t60, t44 * t58, t52 ^ 2 * t56, t44 * t57, t44 ^ 2 / 0.2e1, pkin(1) * t60 + t34 * t44, -pkin(1) * t50 * t66 - t35 * t44 (-t34 * t50 + t35 * t52) * t65, t35 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t56, t61, -t69, -t68, t62, t67, t36, -t15 * t39 + t24 * t30, t16 * t39 + t24 * t32, -t15 * t32 - t16 * t30, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t36, t68, -t67, t61, -t69, t62, t12 * t32 + t13 * t30, -t11 * t30 - t12 * t39, -t11 * t32 + t13 * t39, t11 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t17, -t71, t14, t63, -t70, t26, t10 * t18 + t3 * t29, t10 * t20 - t4 * t29, -t4 * t18 - t3 * t20, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t17, t14, t71, t26, t70, t63, -t1 * t29 + t5 * t18, t1 * t20 - t2 * t18, t2 * t29 - t5 * t20, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
