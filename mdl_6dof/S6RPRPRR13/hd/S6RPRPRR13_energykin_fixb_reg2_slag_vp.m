% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR13_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:25:19
% EndTime: 2019-03-09 04:25:19
% DurationCPUTime: 0.34s
% Computational Cost: add. (1329->79), mult. (4363->177), div. (0->0), fcn. (3551->12), ass. (0->70)
t55 = sin(pkin(12));
t57 = sin(pkin(6));
t58 = cos(pkin(12));
t59 = cos(pkin(6));
t62 = sin(qJ(3));
t80 = cos(pkin(7));
t68 = t62 * t80;
t56 = sin(pkin(7));
t82 = t56 * t62;
t90 = cos(qJ(3));
t36 = ((t90 * t55 + t58 * t68) * t57 + t59 * t82) * qJD(1);
t63 = qJD(1) ^ 2;
t92 = t63 / 0.2e1;
t91 = pkin(3) + pkin(10);
t67 = t80 * t90;
t70 = t56 * t90;
t79 = qJD(1) * t57;
t72 = t58 * t79;
t78 = qJD(1) * t59;
t34 = t62 * t55 * t79 - t67 * t72 - t70 * t78;
t73 = pkin(1) * t78;
t51 = t58 * t73;
t83 = t55 * t57;
t33 = t51 + (pkin(2) * t59 + (-t80 * pkin(9) - qJ(2)) * t83) * qJD(1);
t40 = qJD(2) + (-pkin(9) * t55 * t56 - pkin(2) * t58 - pkin(1)) * t79;
t21 = -t56 * t33 + t80 * t40;
t66 = -t36 * qJ(4) + t21;
t11 = t91 * t34 + t66;
t61 = sin(qJ(5));
t89 = cos(qJ(5));
t43 = t56 * t72 - t80 * t78 - qJD(3);
t69 = qJ(2) * t79;
t46 = t55 * t73 + t58 * t69;
t81 = t57 * t58;
t30 = (t56 * t59 + t80 * t81) * qJD(1) * pkin(9) + t46;
t16 = -t62 * t30 + t33 * t67 + t40 * t70;
t64 = qJD(4) - t16;
t9 = t36 * pkin(4) + t91 * t43 + t64;
t6 = t89 * t11 + t61 * t9;
t88 = cos(qJ(6));
t87 = t36 * t34;
t86 = t36 * t43;
t85 = t43 * t34;
t84 = t57 ^ 2 * t63;
t77 = t34 ^ 2 / 0.2e1;
t76 = t36 ^ 2 / 0.2e1;
t74 = t57 * t59 * t63;
t17 = t90 * t30 + t33 * t68 + t40 * t82;
t71 = t84 / 0.2e1;
t23 = -t89 * t34 - t61 * t43;
t15 = t43 * qJ(4) - t17;
t12 = -t34 * pkin(4) - t15;
t5 = -t61 * t11 + t89 * t9;
t60 = sin(qJ(6));
t52 = -pkin(1) * t79 + qJD(2);
t45 = -t55 * t69 + t51;
t41 = t43 ^ 2 / 0.2e1;
t32 = qJD(5) + t36;
t25 = t61 * t34 - t89 * t43;
t22 = qJD(6) + t23;
t20 = t88 * t25 + t60 * t32;
t18 = t60 * t25 - t88 * t32;
t14 = t43 * pkin(3) + t64;
t13 = t34 * pkin(3) + t66;
t7 = t23 * pkin(5) - t25 * pkin(11) + t12;
t4 = t32 * pkin(11) + t6;
t3 = -t32 * pkin(5) - t5;
t2 = t88 * t4 + t60 * t7;
t1 = -t60 * t4 + t88 * t7;
t8 = [0, 0, 0, 0, 0, t92, 0, 0, 0, 0, t55 ^ 2 * t71, t55 * t58 * t84, t55 * t74, t58 ^ 2 * t71, t58 * t74, t59 ^ 2 * t92 (t45 * t59 - t52 * t81) * qJD(1) (-t46 * t59 + t52 * t83) * qJD(1) (-t45 * t55 + t46 * t58) * t79, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1, t76, -t87, -t86, t77, t85, t41, -t16 * t43 + t21 * t34, t17 * t43 + t21 * t36, -t16 * t36 - t17 * t34, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t41, t86, -t85, t76, -t87, t77, t14 * t36 + t15 * t34, -t13 * t34 - t14 * t43, -t13 * t36 + t15 * t43, t13 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t32, t23 ^ 2 / 0.2e1, -t23 * t32, t32 ^ 2 / 0.2e1, t12 * t23 + t5 * t32, t12 * t25 - t6 * t32, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t22, t18 ^ 2 / 0.2e1, -t18 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t3 * t18, -t2 * t22 + t3 * t20, -t1 * t20 - t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
