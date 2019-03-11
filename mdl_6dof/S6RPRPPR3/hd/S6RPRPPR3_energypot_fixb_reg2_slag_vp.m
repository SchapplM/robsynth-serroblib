% Calculate inertial parameters regressor of potential energy for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:23
% EndTime: 2019-03-09 02:45:23
% DurationCPUTime: 0.13s
% Computational Cost: add. (169->58), mult. (165->68), div. (0->0), fcn. (153->8), ass. (0->33)
t74 = qJ(2) + pkin(6);
t98 = g(3) * t74;
t79 = cos(qJ(3));
t97 = g(3) * t79;
t73 = qJ(1) + pkin(9);
t66 = sin(t73);
t96 = t66 * t79;
t67 = cos(t73);
t95 = t67 * t79;
t75 = sin(qJ(6));
t76 = sin(qJ(3));
t94 = t75 * t76;
t78 = cos(qJ(6));
t93 = t76 * t78;
t92 = qJ(4) * t76;
t91 = t76 * pkin(3) + t74;
t80 = cos(qJ(1));
t90 = t80 * pkin(1) + t67 * pkin(2) + t66 * pkin(7);
t77 = sin(qJ(1));
t89 = t77 * pkin(1) + t66 * pkin(2) - t67 * pkin(7);
t88 = pkin(3) * t95 + t67 * t92 + t90;
t87 = pkin(5) * t76 + pkin(8) * t79;
t86 = g(1) * t67 + g(2) * t66;
t85 = -g(1) * t80 - g(2) * t77;
t84 = -t79 * qJ(4) + t91;
t83 = pkin(3) * t96 + t66 * t92 + t89;
t82 = pkin(4) * t96 + t67 * qJ(5) + t83;
t81 = pkin(4) * t95 - t66 * qJ(5) + t88;
t68 = t76 * pkin(4);
t54 = g(1) * t66 - g(2) * t67;
t53 = g(3) * t76 + t86 * t79;
t52 = t86 * t76 - t97;
t1 = [0, 0, 0, 0, 0, 0, t85, g(1) * t77 - g(2) * t80, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t86, t54, -g(3), t85 * pkin(1) - t98, 0, 0, 0, 0, 0, 0, -t53, t52, -t54, -g(1) * t90 - g(2) * t89 - t98, 0, 0, 0, 0, 0, 0, -t53, -t54, -t52, -g(1) * t88 - g(2) * t83 - g(3) * t84, 0, 0, 0, 0, 0, 0, -t52, t53, t54, -g(1) * t81 - g(2) * t82 - g(3) * (t68 + t84) 0, 0, 0, 0, 0, 0, -g(1) * (-t66 * t75 + t67 * t93) - g(2) * (t66 * t93 + t67 * t75) + t78 * t97, -g(1) * (-t66 * t78 - t67 * t94) - g(2) * (-t66 * t94 + t67 * t78) - t75 * t97, -t53, -g(1) * (t87 * t67 + t81) - g(2) * (t87 * t66 + t82) - g(3) * (t76 * pkin(8) + t68 + (-pkin(5) - qJ(4)) * t79 + t91);];
U_reg  = t1;
