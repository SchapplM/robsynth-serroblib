% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:02
% EndTime: 2019-03-09 03:56:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (149->64), mult. (144->71), div. (0->0), fcn. (132->10), ass. (0->36)
t73 = qJ(3) + pkin(10);
t67 = qJ(5) + t73;
t63 = sin(t67);
t64 = cos(t67);
t98 = pkin(5) * t63 - pkin(9) * t64;
t97 = g(3) * pkin(6);
t96 = pkin(2) + pkin(6);
t93 = g(3) * t64;
t76 = sin(qJ(3));
t92 = t76 * pkin(3);
t75 = sin(qJ(6));
t80 = cos(qJ(1));
t91 = t75 * t80;
t77 = sin(qJ(1));
t90 = t77 * t75;
t78 = cos(qJ(6));
t89 = t77 * t78;
t88 = t78 * t80;
t74 = -qJ(4) - pkin(7);
t87 = t80 * pkin(1) + t77 * qJ(2);
t65 = sin(t73);
t59 = pkin(4) * t65 + t92;
t86 = -qJ(2) - t59;
t79 = cos(qJ(3));
t85 = t79 * pkin(3) + t96;
t84 = t77 * t59 + t87;
t66 = cos(t73);
t83 = pkin(4) * t66 + t85;
t69 = t77 * pkin(1);
t72 = -pkin(8) + t74;
t82 = -t77 * t72 + t69;
t81 = -t80 * qJ(2) + t69;
t60 = g(1) * t77 - g(2) * t80;
t61 = g(1) * t80 + g(2) * t77;
t57 = -g(3) * t63 + t60 * t64;
t1 = [0, 0, 0, 0, 0, 0, -t61, t60, -g(3), -t97, 0, 0, 0, 0, 0, 0, -g(3), t61, -t60, -g(1) * t87 - g(2) * t81 - t97, 0, 0, 0, 0, 0, 0, -g(3) * t79 - t60 * t76, g(3) * t76 - t60 * t79, -t61, -g(1) * (t80 * pkin(7) + t87) - g(2) * (t77 * pkin(7) + t81) - g(3) * t96, 0, 0, 0, 0, 0, 0, -g(3) * t66 - t60 * t65, g(3) * t65 - t60 * t66, -t61, -g(1) * (-t74 * t80 + t77 * t92 + t87) - g(2) * (-t77 * t74 + t69 + (-qJ(2) - t92) * t80) - g(3) * t85, 0, 0, 0, 0, 0, 0, -t60 * t63 - t93, -t57, -t61, -g(1) * (-t80 * t72 + t84) - g(2) * (t86 * t80 + t82) - g(3) * t83, 0, 0, 0, 0, 0, 0, -g(1) * (t63 * t89 + t91) - g(2) * (-t63 * t88 + t90) - t78 * t93, -g(1) * (-t63 * t90 + t88) - g(2) * (t63 * t91 + t89) + t75 * t93, t57, -g(1) * (t98 * t77 + t84) - g(2) * t82 - g(3) * (pkin(5) * t64 + pkin(9) * t63 + t83) + (g(1) * t72 - g(2) * (t86 - t98)) * t80;];
U_reg  = t1;
