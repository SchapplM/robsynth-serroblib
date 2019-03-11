% Calculate inertial parameters regressor of potential energy for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:59:59
% EndTime: 2019-03-09 02:59:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (98->66), mult. (167->71), div. (0->0), fcn. (155->6), ass. (0->34)
t88 = g(3) * pkin(6);
t87 = pkin(2) + pkin(6);
t67 = cos(qJ(3));
t86 = pkin(5) * t67;
t68 = cos(qJ(1));
t85 = g(2) * t68;
t64 = sin(qJ(3));
t84 = g(3) * t64;
t83 = g(1) * qJ(5);
t65 = sin(qJ(1));
t82 = t64 * t65;
t81 = t65 * t67;
t63 = sin(qJ(6));
t80 = t68 * t63;
t66 = cos(qJ(6));
t79 = t68 * t66;
t78 = t68 * pkin(1) + t65 * qJ(2);
t77 = qJ(4) * t67;
t57 = t65 * pkin(7);
t58 = t65 * pkin(1);
t76 = t68 * t77 + t57 + t58;
t75 = t68 * pkin(7) + t78;
t74 = t67 * pkin(3) + t64 * qJ(4) + t87;
t73 = -pkin(3) * t64 - qJ(2);
t72 = -t68 * qJ(2) + t58;
t71 = pkin(3) * t82 + t75;
t70 = t67 * pkin(4) + t74;
t49 = g(1) * t65 - t85;
t69 = -t65 * t77 + t71;
t52 = pkin(4) * t82;
t50 = g(1) * t68 + g(2) * t65;
t48 = t49 * t67 - t84;
t47 = g(1) * t82 + g(3) * t67 - t64 * t85;
t1 = [0, 0, 0, 0, 0, 0, -t50, t49, -g(3), -t88, 0, 0, 0, 0, 0, 0, -g(3), t50, -t49, -g(1) * t78 - g(2) * t72 - t88, 0, 0, 0, 0, 0, 0, -t47, -t48, -t50, -g(1) * t75 - g(2) * (t57 + t72) - g(3) * t87, 0, 0, 0, 0, 0, 0, -t47, -t50, t48, -g(1) * t69 - g(2) * (t73 * t68 + t76) - g(3) * t74, 0, 0, 0, 0, 0, 0, t48, t47, t50, -g(1) * (t52 + t69) - g(2) * (-t65 * qJ(5) + t76) - g(3) * t70 + (t83 - g(2) * (-pkin(4) * t64 + t73)) * t68, 0, 0, 0, 0, 0, 0, -g(1) * (-t66 * t81 - t80) - g(2) * (-t65 * t63 + t67 * t79) - t66 * t84, -g(1) * (t63 * t81 - t79) - g(2) * (-t65 * t66 - t67 * t80) + t63 * t84, -t47, -g(1) * (t52 + t71) - g(2) * t76 - g(3) * (t64 * pkin(5) + t67 * pkin(8) + t70) + (-g(1) * (pkin(8) * t64 - t77 - t86) + g(2) * qJ(5)) * t65 + (t83 + (qJ(2) - t86 - (-pkin(3) - pkin(4) - pkin(8)) * t64) * g(2)) * t68;];
U_reg  = t1;
