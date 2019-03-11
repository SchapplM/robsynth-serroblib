% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_energykin_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:25:02
% EndTime: 2019-03-09 01:25:02
% DurationCPUTime: 0.32s
% Computational Cost: add. (1547->72), mult. (4108->177), div. (0->0), fcn. (3456->16), ass. (0->64)
t60 = cos(pkin(7));
t53 = t60 * qJD(2) + qJD(3);
t56 = sin(pkin(8));
t59 = cos(pkin(8));
t68 = cos(qJ(3));
t57 = sin(pkin(7));
t79 = qJD(2) * t57;
t76 = t68 * t79;
t91 = t53 * t56 + t59 * t76;
t69 = cos(qJ(2));
t58 = sin(pkin(6));
t81 = qJD(1) * t58;
t48 = qJD(2) * pkin(2) + t69 * t81;
t61 = cos(pkin(6));
t80 = qJD(1) * t61;
t90 = t48 * t60 + t57 * t80;
t66 = sin(qJ(2));
t45 = pkin(10) * t79 + t66 * t81;
t65 = sin(qJ(3));
t31 = t68 * t45 + t90 * t65;
t23 = t91 * pkin(11) + t31;
t82 = t90 * t68;
t27 = t53 * pkin(3) + (-pkin(11) * t59 * t79 - t45) * t65 + t82;
t52 = t60 * t80;
t32 = t52 + (-t48 + (-pkin(11) * t56 * t65 - pkin(3) * t68) * qJD(2)) * t57;
t64 = sin(qJ(4));
t67 = cos(qJ(4));
t13 = -t64 * t23 + (t27 * t59 + t32 * t56) * t67;
t83 = t59 * t64;
t84 = t56 * t64;
t14 = t67 * t23 + t27 * t83 + t32 * t84;
t40 = -t59 * t53 + t56 * t76 - qJD(4);
t10 = -t40 * pkin(12) + t14;
t15 = -t56 * t27 + t59 * t32;
t77 = t65 * t79;
t35 = t64 * t77 - t91 * t67;
t37 = t53 * t84 + (t65 * t67 + t68 * t83) * t79;
t12 = t35 * pkin(4) - t37 * pkin(12) + t15;
t63 = sin(qJ(5));
t89 = cos(qJ(5));
t6 = t89 * t10 + t63 * t12;
t88 = cos(qJ(6));
t70 = qJD(2) ^ 2;
t85 = t57 ^ 2 * t70;
t75 = t85 / 0.2e1;
t74 = qJD(2) * t81;
t24 = t63 * t37 + t89 * t40;
t5 = -t63 * t10 + t89 * t12;
t9 = t40 * pkin(4) - t13;
t71 = qJD(1) ^ 2;
t62 = sin(qJ(6));
t39 = -t57 * t48 + t52;
t34 = qJD(5) + t35;
t30 = -t65 * t45 + t82;
t26 = t89 * t37 - t63 * t40;
t22 = qJD(6) + t24;
t18 = t88 * t26 + t62 * t34;
t16 = t62 * t26 - t88 * t34;
t7 = t24 * pkin(5) - t26 * pkin(13) + t9;
t4 = t34 * pkin(13) + t6;
t3 = -t34 * pkin(5) - t5;
t2 = t88 * t4 + t62 * t7;
t1 = -t62 * t4 + t88 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t71 / 0.2e1, 0, 0, 0, 0, 0, t70 / 0.2e1, t69 * t74, -t66 * t74, 0 (t61 ^ 2 / 0.2e1 + (t66 ^ 2 / 0.2e1 + t69 ^ 2 / 0.2e1) * t58 ^ 2) * t71, t65 ^ 2 * t75, t65 * t68 * t85, t53 * t77, t68 ^ 2 * t75, t53 * t76, t53 ^ 2 / 0.2e1, t30 * t53 - t39 * t76, -t31 * t53 + t39 * t77 (-t30 * t65 + t31 * t68) * t79, t31 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, t37 ^ 2 / 0.2e1, -t35 * t37, -t40 * t37, t35 ^ 2 / 0.2e1, t40 * t35, t40 ^ 2 / 0.2e1, -t13 * t40 + t15 * t35, t14 * t40 + t15 * t37, -t13 * t37 - t14 * t35, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, t26 * t34, t24 ^ 2 / 0.2e1, -t24 * t34, t34 ^ 2 / 0.2e1, t9 * t24 + t5 * t34, t9 * t26 - t6 * t34, -t6 * t24 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t22, t16 ^ 2 / 0.2e1, -t16 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t3 * t16, t3 * t18 - t2 * t22, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
