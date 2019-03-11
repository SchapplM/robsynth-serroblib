% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:19:57
% EndTime: 2019-03-09 04:19:57
% DurationCPUTime: 0.14s
% Computational Cost: add. (115->73), mult. (180->88), div. (0->0), fcn. (172->8), ass. (0->39)
t102 = g(3) * pkin(6);
t101 = pkin(2) + pkin(6);
t80 = cos(qJ(1));
t100 = g(2) * t80;
t76 = sin(qJ(3));
t99 = g(3) * t76;
t75 = sin(qJ(5));
t98 = t75 * t76;
t77 = sin(qJ(1));
t97 = t76 * t77;
t79 = cos(qJ(3));
t96 = t77 * t79;
t74 = qJ(5) + qJ(6);
t65 = sin(t74);
t95 = t80 * t65;
t66 = cos(t74);
t94 = t80 * t66;
t93 = t80 * t75;
t78 = cos(qJ(5));
t92 = t80 * t78;
t91 = t80 * pkin(1) + t77 * qJ(2);
t90 = qJ(4) * t79;
t69 = t77 * pkin(7);
t70 = t77 * pkin(1);
t89 = t80 * t90 + t69 + t70;
t88 = t80 * pkin(7) + t91;
t87 = t79 * pkin(3) + t76 * qJ(4) + t101;
t86 = -pkin(3) * t76 - qJ(2);
t85 = -t80 * qJ(2) + t70;
t84 = pkin(3) * t97 + t88;
t59 = g(1) * t77 - t100;
t81 = -pkin(9) - pkin(8);
t83 = -pkin(5) * t75 * t79 - t76 * t81;
t82 = -t77 * t90 + t84;
t64 = t78 * pkin(5) + pkin(4);
t60 = g(1) * t80 + g(2) * t77;
t58 = t59 * t79 - t99;
t57 = g(1) * t97 + g(3) * t79 - t76 * t100;
t1 = [0, 0, 0, 0, 0, 0, -t60, t59, -g(3), -t102, 0, 0, 0, 0, 0, 0, -g(3), t60, -t59, -g(1) * t91 - g(2) * t85 - t102, 0, 0, 0, 0, 0, 0, -t57, -t58, -t60, -g(1) * t88 - g(2) * (t69 + t85) - g(3) * t101, 0, 0, 0, 0, 0, 0, -t60, t57, t58, -g(1) * t82 - g(2) * (t86 * t80 + t89) - g(3) * t87, 0, 0, 0, 0, 0, 0, -g(1) * (-t75 * t96 + t92) - g(2) * (t77 * t78 + t79 * t93) - g(3) * t98, -g(1) * (-t78 * t96 - t93) - g(2) * (-t77 * t75 + t79 * t92) - t78 * t99, -t57, -g(1) * (pkin(8) * t97 + t82) - g(2) * (t77 * pkin(4) + t89) - g(3) * (t79 * pkin(8) + t87) + (-g(1) * pkin(4) - g(2) * (-pkin(8) * t76 + t86)) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (-t65 * t96 + t94) - g(2) * (t77 * t66 + t79 * t95) - t65 * t99, -g(1) * (-t66 * t96 - t95) - g(2) * (-t77 * t65 + t79 * t94) - t66 * t99, -t57, -g(1) * t84 - g(2) * t89 - g(3) * (pkin(5) * t98 - t79 * t81 + t87) + (-g(1) * (t83 - t90) - g(2) * t64) * t77 + (-g(1) * t64 - g(2) * (-t83 + t86)) * t80;];
U_reg  = t1;
