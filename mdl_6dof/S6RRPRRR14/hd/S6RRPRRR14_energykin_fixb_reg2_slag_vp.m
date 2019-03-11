% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR14_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energykin_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:08:36
% EndTime: 2019-03-09 15:08:36
% DurationCPUTime: 0.56s
% Computational Cost: add. (4144->86), mult. (11841->198), div. (0->0), fcn. (10102->16), ass. (0->70)
t91 = cos(pkin(6)) * qJD(1);
t66 = qJD(2) + t91;
t68 = sin(pkin(14));
t73 = cos(pkin(7));
t80 = cos(qJ(2));
t71 = sin(pkin(6));
t92 = qJD(1) * t71;
t87 = t80 * t92;
t84 = t73 * t87;
t70 = sin(pkin(7));
t93 = cos(pkin(14));
t85 = t70 * t93;
t78 = sin(qJ(2));
t88 = t78 * t92;
t47 = -t66 * t85 + t68 * t88 - t93 * t84;
t57 = -t73 * t66 + t70 * t87;
t69 = sin(pkin(8));
t72 = cos(pkin(8));
t105 = t47 * t72 + t57 * t69;
t94 = t68 * t73;
t95 = t68 * t70;
t49 = t66 * t95 + (t93 * t78 + t80 * t94) * t92;
t103 = pkin(11) * t49;
t89 = pkin(1) * t91;
t60 = pkin(10) * t87 + t78 * t89;
t46 = (t66 * t70 + t84) * qJ(3) + t60;
t65 = t80 * t89;
t50 = t66 * pkin(2) + t65 + (-qJ(3) * t73 - pkin(10)) * t88;
t55 = (-qJ(3) * t70 * t78 - pkin(2) * t80 - pkin(1)) * t92;
t31 = t73 * t93 * t50 - t68 * t46 + t55 * t85;
t27 = -t57 * pkin(3) - t72 * t103 + t31;
t38 = -t70 * t50 + t73 * t55 + qJD(3);
t30 = t47 * pkin(3) - t69 * t103 + t38;
t104 = t27 * t72 + t30 * t69;
t32 = t93 * t46 + t50 * t94 + t55 * t95;
t23 = -pkin(11) * t105 + t32;
t77 = sin(qJ(4));
t79 = cos(qJ(4));
t13 = t104 * t79 - t77 * t23;
t14 = t104 * t77 + t79 * t23;
t40 = -t69 * t47 + t72 * t57 - qJD(4);
t10 = -t40 * pkin(12) + t14;
t102 = cos(qJ(5));
t15 = -t69 * t27 + t72 * t30;
t35 = t105 * t79 + t77 * t49;
t37 = -t105 * t77 + t79 * t49;
t12 = t35 * pkin(4) - t37 * pkin(12) + t15;
t76 = sin(qJ(5));
t6 = t102 * t10 + t76 * t12;
t101 = cos(qJ(6));
t81 = qJD(1) ^ 2;
t96 = t71 ^ 2 * t81;
t90 = t80 * t96;
t86 = t96 / 0.2e1;
t24 = t102 * t40 + t76 * t37;
t5 = -t76 * t10 + t102 * t12;
t9 = t40 * pkin(4) - t13;
t75 = sin(qJ(6));
t59 = -pkin(10) * t88 + t65;
t34 = qJD(5) + t35;
t26 = t102 * t37 - t76 * t40;
t22 = qJD(6) + t24;
t18 = t101 * t26 + t75 * t34;
t16 = -t101 * t34 + t75 * t26;
t7 = t24 * pkin(5) - t26 * pkin(13) + t9;
t4 = t34 * pkin(13) + t6;
t3 = -t34 * pkin(5) - t5;
t2 = t101 * t4 + t75 * t7;
t1 = t101 * t7 - t75 * t4;
t8 = [0, 0, 0, 0, 0, t81 / 0.2e1, 0, 0, 0, 0, t78 ^ 2 * t86, t78 * t90, t66 * t88, t80 ^ 2 * t86, t66 * t87, t66 ^ 2 / 0.2e1, pkin(1) * t90 + t59 * t66, -pkin(1) * t78 * t96 - t60 * t66 (-t59 * t78 + t60 * t80) * t92, t60 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t86, t49 ^ 2 / 0.2e1, -t47 * t49, -t57 * t49, t47 ^ 2 / 0.2e1, t47 * t57, t57 ^ 2 / 0.2e1, -t31 * t57 + t38 * t47, t32 * t57 + t38 * t49, -t31 * t49 - t32 * t47, t32 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t37 ^ 2 / 0.2e1, -t35 * t37, -t40 * t37, t35 ^ 2 / 0.2e1, t40 * t35, t40 ^ 2 / 0.2e1, -t13 * t40 + t15 * t35, t14 * t40 + t15 * t37, -t13 * t37 - t14 * t35, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, t26 * t34, t24 ^ 2 / 0.2e1, -t24 * t34, t34 ^ 2 / 0.2e1, t9 * t24 + t5 * t34, t9 * t26 - t6 * t34, -t6 * t24 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t22, t16 ^ 2 / 0.2e1, -t16 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t3 * t16, t3 * t18 - t2 * t22, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
