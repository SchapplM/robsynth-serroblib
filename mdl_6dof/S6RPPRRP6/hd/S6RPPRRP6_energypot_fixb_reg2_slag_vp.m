% Calculate inertial parameters regressor of potential energy for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:17
% EndTime: 2019-03-09 02:11:17
% DurationCPUTime: 0.10s
% Computational Cost: add. (98->51), mult. (164->58), div. (0->0), fcn. (160->6), ass. (0->34)
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t90 = pkin(4) * t65 - pkin(8) * t68;
t89 = g(3) * pkin(6);
t88 = pkin(2) + pkin(6);
t85 = g(3) * t68;
t64 = sin(qJ(5));
t66 = sin(qJ(1));
t84 = t66 * t64;
t67 = cos(qJ(5));
t83 = t66 * t67;
t69 = cos(qJ(1));
t82 = t69 * t64;
t81 = t69 * t67;
t80 = t69 * pkin(1) + t66 * qJ(2);
t79 = pkin(3) + t88;
t78 = t69 * qJ(3) + t80;
t77 = t66 * pkin(1) - t69 * qJ(2);
t76 = t68 * pkin(4) + t65 * pkin(8) + t79;
t75 = t66 * qJ(3) + t77;
t51 = g(1) * t69 + g(2) * t66;
t74 = -t66 * pkin(7) + t78;
t73 = t69 * pkin(7) + t75;
t46 = t65 * t84 - t81;
t48 = t65 * t82 + t83;
t72 = g(1) * t48 + g(2) * t46 + t64 * t85;
t71 = t90 * t69 + t74;
t70 = t90 * t66 + t73;
t50 = g(1) * t66 - g(2) * t69;
t49 = t65 * t81 - t84;
t47 = t65 * t83 + t82;
t45 = -g(3) * t65 + t51 * t68;
t44 = -g(1) * t49 - g(2) * t47 - t67 * t85;
t1 = [0, 0, 0, 0, 0, 0, -t51, t50, -g(3), -t89, 0, 0, 0, 0, 0, 0, -g(3), t51, -t50, -g(1) * t80 - g(2) * t77 - t89, 0, 0, 0, 0, 0, 0, -g(3), -t50, -t51, -g(1) * t78 - g(2) * t75 - g(3) * t88, 0, 0, 0, 0, 0, 0, -t51 * t65 - t85, -t45, t50, -g(1) * t74 - g(2) * t73 - g(3) * t79, 0, 0, 0, 0, 0, 0, t44, t72, t45, -g(1) * t71 - g(2) * t70 - g(3) * t76, 0, 0, 0, 0, 0, 0, t44, t45, -t72, -g(1) * (t49 * pkin(5) + t48 * qJ(6) + t71) - g(2) * (t47 * pkin(5) + t46 * qJ(6) + t70) - g(3) * ((pkin(5) * t67 + qJ(6) * t64) * t68 + t76);];
U_reg  = t1;
