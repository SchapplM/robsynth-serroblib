% Calculate inertial parameters regressor of potential energy for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:58:29
% EndTime: 2019-03-09 06:58:29
% DurationCPUTime: 0.16s
% Computational Cost: add. (211->71), mult. (160->91), div. (0->0), fcn. (152->12), ass. (0->40)
t87 = cos(qJ(5));
t69 = pkin(5) * t87 + pkin(4);
t82 = qJ(3) + qJ(4);
t74 = sin(t82);
t76 = cos(t82);
t90 = -pkin(10) - pkin(9);
t108 = t69 * t76 - t74 * t90;
t107 = g(3) * t74;
t83 = qJ(2) + pkin(6);
t106 = g(3) * t83;
t80 = qJ(1) + pkin(11);
t72 = cos(t80);
t84 = sin(qJ(5));
t104 = t72 * t84;
t81 = qJ(5) + qJ(6);
t73 = sin(t81);
t103 = t73 * t76;
t75 = cos(t81);
t101 = t75 * t76;
t100 = t76 * t84;
t99 = t76 * t87;
t88 = cos(qJ(3));
t70 = pkin(3) * t88 + pkin(2);
t89 = cos(qJ(1));
t79 = t89 * pkin(1);
t98 = t72 * t70 + t79;
t71 = sin(t80);
t86 = sin(qJ(1));
t78 = t86 * pkin(1);
t91 = -pkin(8) - pkin(7);
t97 = t71 * t70 + t72 * t91 + t78;
t85 = sin(qJ(3));
t96 = t85 * pkin(3) + t83;
t95 = -t71 * t91 + t98;
t94 = pkin(4) * t76 + pkin(9) * t74;
t93 = g(1) * t72 + g(2) * t71;
t92 = -g(1) * t89 - g(2) * t86;
t63 = g(1) * t71 - g(2) * t72;
t62 = -g(3) * t76 + t93 * t74;
t1 = [0, 0, 0, 0, 0, 0, t92, g(1) * t86 - g(2) * t89, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t93, t63, -g(3), t92 * pkin(1) - t106, 0, 0, 0, 0, 0, 0, -g(3) * t85 - t93 * t88, -g(3) * t88 + t93 * t85, -t63, -g(1) * (pkin(2) * t72 + pkin(7) * t71 + t79) - g(2) * (pkin(2) * t71 - pkin(7) * t72 + t78) - t106, 0, 0, 0, 0, 0, 0, -t93 * t76 - t107, t62, -t63, -g(1) * t95 - g(2) * t97 - g(3) * t96, 0, 0, 0, 0, 0, 0, -g(1) * (t71 * t84 + t72 * t99) - g(2) * (t71 * t99 - t104) - t87 * t107, -g(1) * (-t72 * t100 + t71 * t87) - g(2) * (-t71 * t100 - t72 * t87) + t84 * t107, -t62, -g(1) * (t94 * t72 + t95) - g(2) * (t94 * t71 + t97) - g(3) * (pkin(4) * t74 - pkin(9) * t76 + t96) 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t101 + t71 * t73) - g(2) * (t71 * t101 - t72 * t73) - t75 * t107, -g(1) * (-t72 * t103 + t71 * t75) - g(2) * (-t71 * t103 - t72 * t75) + t73 * t107, -t62, -g(1) * (t108 * t72 + t98) - g(2) * (-pkin(5) * t104 + t97) - g(3) * (t69 * t74 + t76 * t90 + t96) + (-g(1) * (pkin(5) * t84 - t91) - g(2) * t108) * t71;];
U_reg  = t1;
