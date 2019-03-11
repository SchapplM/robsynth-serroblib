% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:15:10
% EndTime: 2019-03-09 04:15:11
% DurationCPUTime: 0.44s
% Computational Cost: add. (2214->81), mult. (7201->192), div. (0->0), fcn. (6011->14), ass. (0->66)
t63 = sin(pkin(7));
t64 = sin(pkin(6));
t66 = cos(pkin(6));
t65 = cos(pkin(12));
t85 = cos(pkin(7));
t75 = t65 * t85;
t94 = qJD(1) * (t63 * t66 + t64 * t75);
t62 = sin(pkin(12));
t83 = qJD(1) * t64;
t76 = qJ(2) * t83;
t82 = qJD(1) * t66;
t78 = pkin(1) * t82;
t52 = t62 * t78 + t65 * t76;
t39 = pkin(9) * t94 + t52;
t47 = qJD(2) + (-pkin(9) * t62 * t63 - pkin(2) * t65 - pkin(1)) * t83;
t69 = sin(qJ(3));
t70 = cos(qJ(3));
t57 = t65 * t78;
t88 = t62 * t64;
t41 = t57 + (pkin(2) * t66 + (-t85 * pkin(9) - qJ(2)) * t88) * qJD(1);
t74 = t85 * t41;
t28 = -t69 * t39 + (t47 * t63 + t74) * t70;
t71 = qJD(1) ^ 2;
t92 = t71 / 0.2e1;
t30 = -t63 * t41 + t85 * t47;
t42 = t69 * t62 * t83 - t70 * t94;
t86 = t63 * t69;
t44 = (t66 * t86 + (t62 * t70 + t69 * t75) * t64) * qJD(1);
t21 = t42 * pkin(3) - t44 * qJ(4) + t30;
t29 = t70 * t39 + t47 * t86 + t69 * t74;
t49 = t63 * t65 * t83 - t85 * t82 - qJD(3);
t23 = -t49 * qJ(4) + t29;
t61 = sin(pkin(13));
t84 = cos(pkin(13));
t13 = t61 * t21 + t84 * t23;
t32 = t61 * t44 + t84 * t49;
t11 = -t32 * pkin(10) + t13;
t68 = sin(qJ(5));
t12 = t84 * t21 - t61 * t23;
t34 = t84 * t44 - t61 * t49;
t9 = t42 * pkin(4) - t34 * pkin(10) + t12;
t91 = cos(qJ(5));
t6 = t91 * t11 + t68 * t9;
t90 = cos(qJ(6));
t89 = t64 ^ 2 * t71;
t80 = t42 ^ 2 / 0.2e1;
t79 = t64 * t66 * t71;
t77 = t89 / 0.2e1;
t25 = t91 * t32 + t68 * t34;
t5 = -t68 * t11 + t91 * t9;
t22 = t49 * pkin(3) + qJD(4) - t28;
t14 = t32 * pkin(4) + t22;
t67 = sin(qJ(6));
t58 = -pkin(1) * t83 + qJD(2);
t51 = -t62 * t76 + t57;
t40 = qJD(5) + t42;
t27 = -t68 * t32 + t91 * t34;
t24 = qJD(6) + t25;
t17 = t90 * t27 + t67 * t40;
t15 = t67 * t27 - t90 * t40;
t7 = t25 * pkin(5) - t27 * pkin(11) + t14;
t4 = t40 * pkin(11) + t6;
t3 = -t40 * pkin(5) - t5;
t2 = t90 * t4 + t67 * t7;
t1 = -t67 * t4 + t90 * t7;
t8 = [0, 0, 0, 0, 0, t92, 0, 0, 0, 0, t62 ^ 2 * t77, t62 * t65 * t89, t62 * t79, t65 ^ 2 * t77, t65 * t79, t66 ^ 2 * t92 (-t58 * t64 * t65 + t51 * t66) * qJD(1) (-t52 * t66 + t58 * t88) * qJD(1) (-t51 * t62 + t52 * t65) * t83, t52 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1 + t58 ^ 2 / 0.2e1, t44 ^ 2 / 0.2e1, -t44 * t42, -t44 * t49, t80, t42 * t49, t49 ^ 2 / 0.2e1, -t28 * t49 + t30 * t42, t29 * t49 + t30 * t44, -t28 * t44 - t29 * t42, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t42, t32 ^ 2 / 0.2e1, -t32 * t42, t80, t12 * t42 + t22 * t32, -t13 * t42 + t22 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t40, t25 ^ 2 / 0.2e1, -t25 * t40, t40 ^ 2 / 0.2e1, t14 * t25 + t5 * t40, t14 * t27 - t6 * t40, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t24, t15 ^ 2 / 0.2e1, -t15 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t3 * t15, t3 * t17 - t2 * t24, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
