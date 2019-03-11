% Calculate inertial parameters regressor of potential energy for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:47
% EndTime: 2019-03-09 02:23:48
% DurationCPUTime: 0.17s
% Computational Cost: add. (171->66), mult. (149->83), div. (0->0), fcn. (141->10), ass. (0->37)
t74 = cos(qJ(5));
t61 = t74 * pkin(5) + pkin(4);
t72 = sin(qJ(4));
t75 = cos(qJ(4));
t77 = -pkin(9) - pkin(8);
t96 = t61 * t72 + t75 * t77;
t68 = qJ(1) + pkin(10);
t63 = cos(t68);
t95 = g(2) * t63;
t70 = qJ(2) + pkin(6);
t94 = g(3) * t70;
t93 = g(3) * t75;
t62 = sin(t68);
t71 = sin(qJ(5));
t91 = t62 * t71;
t69 = qJ(5) + qJ(6);
t64 = sin(t69);
t90 = t64 * t72;
t65 = cos(t69);
t89 = t65 * t72;
t88 = t71 * t72;
t87 = t72 * t74;
t73 = sin(qJ(1));
t85 = t73 * pkin(1) + t62 * pkin(2);
t84 = pkin(3) + t70;
t57 = t62 * pkin(7);
t83 = t57 + t85;
t76 = cos(qJ(1));
t82 = t76 * pkin(1) + t63 * pkin(2) + t62 * qJ(3);
t81 = t63 * pkin(7) + t82;
t80 = -t63 * qJ(3) + t85;
t79 = pkin(4) * t72 - pkin(8) * t75;
t54 = g(1) * t62 - t95;
t78 = -g(1) * t76 - g(2) * t73;
t55 = g(1) * t63 + g(2) * t62;
t53 = -g(3) * t72 + t54 * t75;
t1 = [0, 0, 0, 0, 0, 0, t78, g(1) * t73 - g(2) * t76, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t55, t54, -g(3), t78 * pkin(1) - t94, 0, 0, 0, 0, 0, 0, -g(3), t55, -t54, -g(1) * t82 - g(2) * t80 - t94, 0, 0, 0, 0, 0, 0, -t54 * t72 - t93, -t53, -t55, -g(1) * t81 - g(2) * (t57 + t80) - g(3) * t84, 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t87 + t63 * t71) - g(2) * (-t63 * t87 + t91) - t74 * t93, -g(1) * (-t62 * t88 + t63 * t74) - g(2) * (t62 * t74 + t63 * t88) + t71 * t93, t53, -g(1) * (t79 * t62 + t81) - g(2) * t83 - g(3) * (t75 * pkin(4) + t72 * pkin(8) + t84) - (-qJ(3) - t79) * t95, 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t89 + t63 * t64) - g(2) * (t62 * t64 - t63 * t89) - t65 * t93, -g(1) * (-t62 * t90 + t63 * t65) - g(2) * (t62 * t65 + t63 * t90) + t64 * t93, t53, -g(1) * (t96 * t62 + t81) - g(2) * (pkin(5) * t91 + t83) - g(3) * (t75 * t61 - t72 * t77 + t84) + (-g(1) * pkin(5) * t71 - g(2) * (-qJ(3) - t96)) * t63;];
U_reg  = t1;
