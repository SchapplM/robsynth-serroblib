% Calculate minimal parameter regressor of potential energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:25:59
% EndTime: 2022-01-23 09:25:59
% DurationCPUTime: 0.09s
% Computational Cost: add. (83->46), mult. (104->70), div. (0->0), fcn. (110->10), ass. (0->32)
t76 = sin(pkin(8));
t97 = g(3) * t76;
t75 = qJ(3) + pkin(9);
t74 = qJ(5) + t75;
t69 = sin(t74);
t80 = sin(qJ(1));
t96 = t80 * t69;
t70 = cos(t74);
t95 = t80 * t70;
t72 = sin(t75);
t94 = t80 * t72;
t73 = cos(t75);
t93 = t80 * t73;
t79 = sin(qJ(3));
t92 = t80 * t79;
t81 = cos(qJ(3));
t91 = t80 * t81;
t82 = cos(qJ(1));
t90 = t82 * t69;
t89 = t82 * t70;
t88 = t82 * t72;
t87 = t82 * t73;
t86 = t82 * t79;
t85 = t82 * t81;
t84 = pkin(3) * t81 + pkin(2);
t83 = -g(1) * t82 - g(2) * t80;
t78 = qJ(4) + pkin(6);
t77 = cos(pkin(8));
t71 = t79 * pkin(3) + qJ(2);
t68 = g(1) * t80 - g(2) * t82;
t67 = t78 * t76 + t84 * t77 + pkin(1);
t1 = [0, t83, t68, t83 * t77 - t97, -t68, -g(1) * (t82 * pkin(1) + t80 * qJ(2)) - g(2) * (t80 * pkin(1) - t82 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(1) * (t77 * t85 + t92) - g(2) * (t77 * t91 - t86) - t81 * t97, -g(1) * (-t77 * t86 + t91) - g(2) * (-t77 * t92 - t85) + t79 * t97, -g(1) * (t77 * t87 + t94) - g(2) * (t77 * t93 - t88) - t73 * t97, -g(1) * (-t77 * t88 + t93) - g(2) * (-t77 * t94 - t87) + t72 * t97, g(3) * t77 + t83 * t76, -g(1) * (t67 * t82 + t71 * t80) - g(2) * (t67 * t80 - t71 * t82) - g(3) * (t84 * t76 - t77 * t78 + pkin(5)), 0, 0, 0, 0, 0, -g(1) * (t77 * t89 + t96) - g(2) * (t77 * t95 - t90) - t70 * t97, -g(1) * (-t77 * t90 + t95) - g(2) * (-t77 * t96 - t89) + t69 * t97;];
U_reg = t1;
