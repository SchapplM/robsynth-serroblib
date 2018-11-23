% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energykin_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:16:21
% EndTime: 2018-11-23 11:16:22
% DurationCPUTime: 0.61s
% Computational Cost: add. (4343->86), mult. (11841->198), div. (0->0), fcn. (10102->16), ass. (0->70)
t102 = cos(qJ(3));
t91 = cos(pkin(6)) * qJD(1);
t66 = qJD(2) + t91;
t77 = sin(qJ(3));
t72 = cos(pkin(7));
t80 = cos(qJ(2));
t70 = sin(pkin(6));
t92 = qJD(1) * t70;
t87 = t80 * t92;
t84 = t72 * t87;
t69 = sin(pkin(7));
t85 = t69 * t102;
t78 = sin(qJ(2));
t88 = t78 * t92;
t48 = -t102 * t84 - t66 * t85 + t77 * t88;
t57 = -t72 * t66 + t69 * t87 - qJD(3);
t68 = sin(pkin(8));
t71 = cos(pkin(8));
t105 = t48 * t71 + t57 * t68;
t93 = t72 * t77;
t94 = t69 * t77;
t50 = t66 * t94 + (t102 * t78 + t80 * t93) * t92;
t103 = pkin(12) * t50;
t89 = pkin(1) * t91;
t60 = pkin(10) * t87 + t78 * t89;
t46 = (t66 * t69 + t84) * pkin(11) + t60;
t65 = t80 * t89;
t47 = t66 * pkin(2) + t65 + (-pkin(11) * t72 - pkin(10)) * t88;
t56 = (-pkin(11) * t69 * t78 - pkin(2) * t80 - pkin(1)) * t92;
t31 = t72 * t102 * t47 - t77 * t46 + t56 * t85;
t27 = -t57 * pkin(3) - t71 * t103 + t31;
t39 = -t69 * t47 + t72 * t56;
t30 = t48 * pkin(3) - t68 * t103 + t39;
t104 = t27 * t71 + t30 * t68;
t32 = t102 * t46 + t47 * t93 + t56 * t94;
t26 = -pkin(12) * t105 + t32;
t76 = sin(qJ(4));
t79 = cos(qJ(4));
t13 = t104 * t79 - t76 * t26;
t14 = t104 * t76 + t79 * t26;
t40 = -t68 * t48 + t71 * t57 - qJD(4);
t10 = -t40 * pkin(13) + t14;
t101 = cos(qJ(5));
t15 = -t68 * t27 + t71 * t30;
t35 = t105 * t79 + t76 * t50;
t37 = -t105 * t76 + t79 * t50;
t12 = t35 * pkin(4) - t37 * pkin(13) + t15;
t75 = sin(qJ(5));
t6 = t101 * t10 + t75 * t12;
t100 = cos(qJ(6));
t81 = qJD(1) ^ 2;
t95 = t70 ^ 2 * t81;
t90 = t80 * t95;
t86 = t95 / 0.2e1;
t23 = t101 * t40 + t75 * t37;
t5 = -t75 * t10 + t101 * t12;
t9 = t40 * pkin(4) - t13;
t74 = sin(qJ(6));
t59 = -pkin(10) * t88 + t65;
t34 = qJD(5) + t35;
t25 = t101 * t37 - t75 * t40;
t22 = qJD(6) + t23;
t18 = t100 * t25 + t74 * t34;
t16 = -t100 * t34 + t74 * t25;
t7 = t23 * pkin(5) - t25 * pkin(14) + t9;
t4 = t34 * pkin(14) + t6;
t3 = -t34 * pkin(5) - t5;
t2 = t100 * t4 + t74 * t7;
t1 = t100 * t7 - t74 * t4;
t8 = [0, 0, 0, 0, 0, t81 / 0.2e1, 0, 0, 0, 0, t78 ^ 2 * t86, t78 * t90, t66 * t88, t80 ^ 2 * t86, t66 * t87, t66 ^ 2 / 0.2e1, pkin(1) * t90 + t59 * t66, -pkin(1) * t78 * t95 - t60 * t66 (-t59 * t78 + t60 * t80) * t92, t60 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t86, t50 ^ 2 / 0.2e1, -t48 * t50, -t57 * t50, t48 ^ 2 / 0.2e1, t48 * t57, t57 ^ 2 / 0.2e1, -t31 * t57 + t39 * t48, t32 * t57 + t39 * t50, -t31 * t50 - t32 * t48, t32 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, t37 ^ 2 / 0.2e1, -t35 * t37, -t40 * t37, t35 ^ 2 / 0.2e1, t40 * t35, t40 ^ 2 / 0.2e1, -t13 * t40 + t15 * t35, t14 * t40 + t15 * t37, -t13 * t37 - t14 * t35, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t34, t23 ^ 2 / 0.2e1, -t23 * t34, t34 ^ 2 / 0.2e1, t9 * t23 + t5 * t34, t9 * t25 - t6 * t34, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t22, t16 ^ 2 / 0.2e1, -t16 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t3 * t16, t3 * t18 - t2 * t22, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
