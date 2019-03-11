% Calculate inertial parameters regressor of potential energy for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:17:54
% EndTime: 2019-03-09 07:17:54
% DurationCPUTime: 0.16s
% Computational Cost: add. (149->64), mult. (144->71), div. (0->0), fcn. (132->10), ass. (0->36)
t73 = qJ(3) + qJ(4);
t68 = qJ(5) + t73;
t63 = sin(t68);
t64 = cos(t68);
t98 = pkin(5) * t63 - pkin(10) * t64;
t97 = g(3) * pkin(6);
t96 = pkin(2) + pkin(6);
t80 = -pkin(8) - pkin(7);
t93 = g(3) * t64;
t75 = sin(qJ(3));
t92 = t75 * pkin(3);
t74 = sin(qJ(6));
t76 = sin(qJ(1));
t91 = t76 * t74;
t77 = cos(qJ(6));
t90 = t76 * t77;
t79 = cos(qJ(1));
t89 = t79 * t74;
t88 = t79 * t77;
t87 = t79 * pkin(1) + t76 * qJ(2);
t65 = sin(t73);
t59 = pkin(4) * t65 + t92;
t86 = -qJ(2) - t59;
t78 = cos(qJ(3));
t85 = t78 * pkin(3) + t96;
t84 = t76 * t59 + t87;
t66 = cos(t73);
t83 = pkin(4) * t66 + t85;
t69 = t76 * pkin(1);
t72 = -pkin(9) + t80;
t82 = -t76 * t72 + t69;
t81 = -t79 * qJ(2) + t69;
t60 = g(1) * t76 - g(2) * t79;
t61 = g(1) * t79 + g(2) * t76;
t57 = -g(3) * t63 + t60 * t64;
t1 = [0, 0, 0, 0, 0, 0, -t61, t60, -g(3), -t97, 0, 0, 0, 0, 0, 0, -g(3), t61, -t60, -g(1) * t87 - g(2) * t81 - t97, 0, 0, 0, 0, 0, 0, -g(3) * t78 - t60 * t75, g(3) * t75 - t60 * t78, -t61, -g(1) * (t79 * pkin(7) + t87) - g(2) * (t76 * pkin(7) + t81) - g(3) * t96, 0, 0, 0, 0, 0, 0, -g(3) * t66 - t60 * t65, g(3) * t65 - t60 * t66, -t61, -g(1) * (t76 * t92 - t79 * t80 + t87) - g(2) * (-t76 * t80 + t69 + (-qJ(2) - t92) * t79) - g(3) * t85, 0, 0, 0, 0, 0, 0, -t60 * t63 - t93, -t57, -t61, -g(1) * (-t79 * t72 + t84) - g(2) * (t79 * t86 + t82) - g(3) * t83, 0, 0, 0, 0, 0, 0, -g(1) * (t63 * t90 + t89) - g(2) * (-t63 * t88 + t91) - t77 * t93, -g(1) * (-t63 * t91 + t88) - g(2) * (t63 * t89 + t90) + t74 * t93, t57, -g(1) * (t98 * t76 + t84) - g(2) * t82 - g(3) * (t64 * pkin(5) + t63 * pkin(10) + t83) + (g(1) * t72 - g(2) * (t86 - t98)) * t79;];
U_reg  = t1;
