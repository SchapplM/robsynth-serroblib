% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:10:15
% EndTime: 2019-03-09 19:10:15
% DurationCPUTime: 0.34s
% Computational Cost: add. (2577->81), mult. (7223->185), div. (0->0), fcn. (6036->14), ass. (0->65)
t71 = sin(qJ(2));
t72 = cos(qJ(2));
t65 = sin(pkin(6));
t84 = qJD(1) * t65;
t79 = t72 * t84;
t83 = cos(pkin(6)) * qJD(1);
t81 = pkin(1) * t83;
t54 = pkin(9) * t79 + t71 * t81;
t61 = qJD(2) + t83;
t64 = sin(pkin(7));
t85 = cos(pkin(7));
t75 = t72 * t85;
t74 = t75 * t84;
t40 = (t61 * t64 + t74) * pkin(10) + t54;
t49 = (-pkin(10) * t64 * t71 - pkin(2) * t72 - pkin(1)) * t84;
t70 = sin(qJ(3));
t60 = t72 * t81;
t80 = t71 * t84;
t41 = t61 * pkin(2) + t60 + (-t85 * pkin(10) - pkin(9)) * t80;
t76 = t41 * t85;
t90 = cos(qJ(3));
t77 = t64 * t90;
t23 = -t70 * t40 + t49 * t77 + t90 * t76;
t86 = t64 * t70;
t44 = t61 * t86 + (t70 * t75 + t90 * t71) * t84;
t51 = -t85 * t61 + t64 * t79 - qJD(3);
t19 = -t51 * pkin(3) - t44 * qJ(4) + t23;
t24 = t90 * t40 + t49 * t86 + t70 * t76;
t42 = -t61 * t77 + t70 * t80 - t90 * t74;
t22 = -t42 * qJ(4) + t24;
t63 = sin(pkin(13));
t66 = cos(pkin(13));
t12 = t63 * t19 + t66 * t22;
t10 = -t51 * pkin(11) + t12;
t35 = -t64 * t41 + t85 * t49;
t29 = t42 * pkin(3) + qJD(4) + t35;
t32 = t66 * t42 + t63 * t44;
t34 = -t63 * t42 + t66 * t44;
t14 = t32 * pkin(4) - t34 * pkin(11) + t29;
t69 = sin(qJ(5));
t89 = cos(qJ(5));
t6 = t89 * t10 + t69 * t14;
t88 = cos(qJ(6));
t73 = qJD(1) ^ 2;
t87 = t65 ^ 2 * t73;
t82 = t72 * t87;
t78 = t87 / 0.2e1;
t11 = t66 * t19 - t63 * t22;
t26 = t69 * t34 + t89 * t51;
t9 = t51 * pkin(4) - t11;
t5 = -t69 * t10 + t89 * t14;
t68 = sin(qJ(6));
t53 = -pkin(9) * t80 + t60;
t50 = t51 ^ 2 / 0.2e1;
t31 = qJD(5) + t32;
t28 = t89 * t34 - t69 * t51;
t25 = qJD(6) + t26;
t17 = t88 * t28 + t68 * t31;
t15 = t68 * t28 - t88 * t31;
t7 = t26 * pkin(5) - t28 * pkin(12) + t9;
t4 = t31 * pkin(12) + t6;
t3 = -t31 * pkin(5) - t5;
t2 = t88 * t4 + t68 * t7;
t1 = -t68 * t4 + t88 * t7;
t8 = [0, 0, 0, 0, 0, t73 / 0.2e1, 0, 0, 0, 0, t71 ^ 2 * t78, t71 * t82, t61 * t80, t72 ^ 2 * t78, t61 * t79, t61 ^ 2 / 0.2e1, pkin(1) * t82 + t53 * t61, -pkin(1) * t71 * t87 - t54 * t61 (-t53 * t71 + t54 * t72) * t84, t54 ^ 2 / 0.2e1 + t53 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t78, t44 ^ 2 / 0.2e1, -t44 * t42, -t44 * t51, t42 ^ 2 / 0.2e1, t42 * t51, t50, -t23 * t51 + t35 * t42, t24 * t51 + t35 * t44, -t23 * t44 - t24 * t42, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, -t34 * t51, t32 ^ 2 / 0.2e1, t32 * t51, t50, -t11 * t51 + t29 * t32, t12 * t51 + t29 * t34, -t11 * t34 - t12 * t32, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * t31, t26 ^ 2 / 0.2e1, -t26 * t31, t31 ^ 2 / 0.2e1, t9 * t26 + t5 * t31, t9 * t28 - t6 * t31, -t6 * t26 - t5 * t28, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t25, t15 ^ 2 / 0.2e1, -t15 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t3 * t15, t3 * t17 - t2 * t25, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
