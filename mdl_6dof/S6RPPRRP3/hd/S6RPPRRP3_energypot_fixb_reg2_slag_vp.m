% Calculate inertial parameters regressor of potential energy for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:47
% EndTime: 2019-03-09 02:03:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (169->54), mult. (162->63), div. (0->0), fcn. (158->8), ass. (0->35)
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t95 = pkin(4) * t74 - pkin(8) * t77;
t72 = qJ(2) + pkin(6);
t92 = g(3) * t72;
t91 = g(3) * t77;
t73 = sin(qJ(5));
t90 = t73 * t74;
t76 = cos(qJ(5));
t89 = t74 * t76;
t71 = qJ(1) + pkin(9);
t65 = sin(t71);
t75 = sin(qJ(1));
t88 = t75 * pkin(1) + t65 * pkin(2);
t87 = pkin(3) + t72;
t66 = cos(t71);
t78 = cos(qJ(1));
t86 = t78 * pkin(1) + t66 * pkin(2) + t65 * qJ(3);
t85 = t66 * pkin(7) + t86;
t84 = t77 * pkin(4) + t74 * pkin(8) + t87;
t83 = -t66 * qJ(3) + t88;
t55 = g(1) * t65 - g(2) * t66;
t82 = -g(1) * t78 - g(2) * t75;
t51 = t65 * t90 - t66 * t76;
t53 = t65 * t76 + t66 * t90;
t81 = g(1) * t51 - g(2) * t53 + t73 * t91;
t80 = t95 * t65 + t85;
t61 = t65 * pkin(7);
t79 = t61 + t88 + (-qJ(3) - t95) * t66;
t56 = g(1) * t66 + g(2) * t65;
t54 = t65 * t73 - t66 * t89;
t52 = t65 * t89 + t66 * t73;
t50 = -g(3) * t74 + t55 * t77;
t49 = -g(1) * t52 - g(2) * t54 - t76 * t91;
t1 = [0, 0, 0, 0, 0, 0, t82, g(1) * t75 - g(2) * t78, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t56, t55, -g(3), t82 * pkin(1) - t92, 0, 0, 0, 0, 0, 0, -g(3), t56, -t55, -g(1) * t86 - g(2) * t83 - t92, 0, 0, 0, 0, 0, 0, -t55 * t74 - t91, -t50, -t56, -g(1) * t85 - g(2) * (t61 + t83) - g(3) * t87, 0, 0, 0, 0, 0, 0, t49, t81, t50, -g(1) * t80 - g(2) * t79 - g(3) * t84, 0, 0, 0, 0, 0, 0, t49, t50, -t81, -g(1) * (t52 * pkin(5) + t51 * qJ(6) + t80) - g(2) * (t54 * pkin(5) - t53 * qJ(6) + t79) - g(3) * ((pkin(5) * t76 + qJ(6) * t73) * t77 + t84);];
U_reg  = t1;
