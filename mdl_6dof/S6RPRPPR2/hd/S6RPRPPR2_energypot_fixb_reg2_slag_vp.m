% Calculate inertial parameters regressor of potential energy for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:38
% EndTime: 2019-03-09 02:42:38
% DurationCPUTime: 0.15s
% Computational Cost: add. (191->62), mult. (147->70), div. (0->0), fcn. (135->10), ass. (0->39)
t81 = qJ(3) + pkin(10);
t76 = cos(t81);
t82 = qJ(1) + pkin(9);
t77 = cos(t82);
t102 = t76 * t77;
t74 = sin(t81);
t98 = qJ(5) * t74;
t108 = pkin(4) * t102 + t77 * t98;
t107 = g(3) * t76;
t84 = qJ(2) + pkin(6);
t106 = g(3) * t84;
t75 = sin(t82);
t105 = t75 * t76;
t85 = sin(qJ(6));
t104 = t75 * t85;
t88 = cos(qJ(6));
t103 = t75 * t88;
t101 = t77 * t85;
t100 = t77 * t88;
t89 = cos(qJ(3));
t73 = pkin(3) * t89 + pkin(2);
t90 = cos(qJ(1));
t80 = t90 * pkin(1);
t99 = t77 * t73 + t80;
t87 = sin(qJ(1));
t79 = t87 * pkin(1);
t83 = -qJ(4) - pkin(7);
t97 = t75 * t73 + t77 * t83 + t79;
t86 = sin(qJ(3));
t96 = t86 * pkin(3) + t84;
t95 = -t75 * t83 + t99;
t94 = pkin(4) * t105 + t75 * t98 + t97;
t93 = g(1) * t77 + g(2) * t75;
t92 = -g(1) * t90 - g(2) * t87;
t91 = t74 * pkin(4) - qJ(5) * t76 + t96;
t63 = g(1) * t75 - g(2) * t77;
t62 = g(3) * t74 + t93 * t76;
t61 = t93 * t74 - t107;
t1 = [0, 0, 0, 0, 0, 0, t92, g(1) * t87 - g(2) * t90, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t93, t63, -g(3), t92 * pkin(1) - t106, 0, 0, 0, 0, 0, 0, -g(3) * t86 - t93 * t89, -g(3) * t89 + t93 * t86, -t63, -g(1) * (pkin(2) * t77 + pkin(7) * t75 + t80) - g(2) * (pkin(2) * t75 - pkin(7) * t77 + t79) - t106, 0, 0, 0, 0, 0, 0, -t62, t61, -t63, -g(1) * t95 - g(2) * t97 - g(3) * t96, 0, 0, 0, 0, 0, 0, -t63, t62, -t61, -g(1) * (t95 + t108) - g(2) * t94 - g(3) * t91, 0, 0, 0, 0, 0, 0, -g(1) * (t74 * t101 + t103) - g(2) * (t74 * t104 - t100) + t85 * t107, -g(1) * (t74 * t100 - t104) - g(2) * (t74 * t103 + t101) + t88 * t107, -t62, -g(1) * (pkin(8) * t102 + (pkin(5) - t83) * t75 + t99 + t108) - g(2) * (-pkin(5) * t77 + pkin(8) * t105 + t94) - g(3) * (pkin(8) * t74 + t91);];
U_reg  = t1;
