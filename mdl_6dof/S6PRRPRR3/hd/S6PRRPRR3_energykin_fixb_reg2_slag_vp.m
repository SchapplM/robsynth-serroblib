% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:07:41
% EndTime: 2019-03-08 22:07:42
% DurationCPUTime: 0.25s
% Computational Cost: add. (859->65), mult. (2360->159), div. (0->0), fcn. (1918->14), ass. (0->60)
t61 = cos(qJ(2));
t52 = sin(pkin(6));
t73 = qJD(1) * t52;
t40 = qJD(2) * pkin(2) + t61 * t73;
t51 = sin(pkin(7));
t54 = cos(pkin(7));
t55 = cos(pkin(6));
t72 = qJD(1) * t55;
t79 = t40 * t54 + t51 * t72;
t59 = sin(qJ(2));
t71 = qJD(2) * t51;
t39 = pkin(9) * t71 + t59 * t73;
t47 = t54 * qJD(2) + qJD(3);
t58 = sin(qJ(3));
t65 = qJ(4) * t71;
t60 = cos(qJ(3));
t74 = t79 * t60;
t19 = t47 * pkin(3) + (-t39 - t65) * t58 + t74;
t24 = t60 * t39 + t79 * t58;
t22 = t60 * t65 + t24;
t50 = sin(pkin(13));
t53 = cos(pkin(13));
t12 = t50 * t19 + t53 * t22;
t10 = t47 * pkin(10) + t12;
t46 = t54 * t72;
t70 = qJD(2) * t60;
t29 = qJD(4) + t46 + (-pkin(3) * t70 - t40) * t51;
t67 = t51 * t70;
t68 = t58 * t71;
t33 = t50 * t68 - t53 * t67;
t35 = (t50 * t60 + t53 * t58) * t71;
t14 = t33 * pkin(4) - t35 * pkin(10) + t29;
t57 = sin(qJ(5));
t78 = cos(qJ(5));
t6 = t78 * t10 + t57 * t14;
t77 = cos(qJ(6));
t62 = qJD(2) ^ 2;
t75 = t51 ^ 2 * t62;
t66 = t75 / 0.2e1;
t64 = qJD(2) * t73;
t11 = t53 * t19 - t50 * t22;
t26 = t57 * t35 - t78 * t47;
t9 = -t47 * pkin(4) - t11;
t5 = -t57 * t10 + t78 * t14;
t63 = qJD(1) ^ 2;
t56 = sin(qJ(6));
t45 = t47 ^ 2 / 0.2e1;
t32 = qJD(5) + t33;
t31 = -t51 * t40 + t46;
t28 = t78 * t35 + t57 * t47;
t25 = qJD(6) + t26;
t23 = -t58 * t39 + t74;
t18 = t77 * t28 + t56 * t32;
t16 = t56 * t28 - t77 * t32;
t7 = t26 * pkin(5) - t28 * pkin(11) + t9;
t4 = t32 * pkin(11) + t6;
t3 = -t32 * pkin(5) - t5;
t2 = t77 * t4 + t56 * t7;
t1 = -t56 * t4 + t77 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t63 / 0.2e1, 0, 0, 0, 0, 0, t62 / 0.2e1, t61 * t64, -t59 * t64, 0 (t55 ^ 2 / 0.2e1 + (t59 ^ 2 / 0.2e1 + t61 ^ 2 / 0.2e1) * t52 ^ 2) * t63, t58 ^ 2 * t66, t58 * t60 * t75, t47 * t68, t60 ^ 2 * t66, t47 * t67, t45, t23 * t47 - t31 * t67, -t24 * t47 + t31 * t68 (-t23 * t58 + t24 * t60) * t71, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * t47, t33 ^ 2 / 0.2e1, -t33 * t47, t45, t11 * t47 + t29 * t33, -t12 * t47 + t29 * t35, -t11 * t35 - t12 * t33, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * t32, t26 ^ 2 / 0.2e1, -t26 * t32, t32 ^ 2 / 0.2e1, t9 * t26 + t5 * t32, t9 * t28 - t6 * t32, -t6 * t26 - t5 * t28, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t25, t16 ^ 2 / 0.2e1, -t16 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t3 * t16, t3 * t18 - t2 * t25, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
