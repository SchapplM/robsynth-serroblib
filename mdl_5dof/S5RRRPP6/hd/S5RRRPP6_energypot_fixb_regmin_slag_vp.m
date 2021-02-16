% Calculate minimal parameter regressor of potential energy for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:38
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:37:05
% EndTime: 2021-01-15 22:37:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (110->49), mult. (145->67), div. (0->0), fcn. (151->10), ass. (0->33)
t84 = sin(qJ(2));
t99 = g(3) * t84;
t82 = qJ(4) + pkin(7);
t98 = t82 * t84 + pkin(1);
t87 = cos(qJ(2));
t97 = -t87 * t82 + pkin(5);
t85 = sin(qJ(1));
t96 = t85 * t87;
t79 = qJ(3) + pkin(8);
t76 = sin(t79);
t88 = cos(qJ(1));
t95 = t88 * t76;
t77 = cos(t79);
t94 = t88 * t77;
t83 = sin(qJ(3));
t93 = t88 * t83;
t86 = cos(qJ(3));
t92 = t88 * t86;
t91 = g(1) * t88 + g(2) * t85;
t80 = sin(pkin(8));
t81 = cos(pkin(8));
t70 = pkin(4) * t81 + qJ(5) * t80 + pkin(3);
t71 = -t80 * pkin(4) + qJ(5) * t81;
t90 = t70 * t83 - t71 * t86 + pkin(6);
t89 = g(1) * (-t85 * t77 + t87 * t95) + g(2) * (t76 * t96 + t94) + t76 * t99;
t75 = t86 * pkin(3) + pkin(2);
t74 = t83 * pkin(3) + pkin(6);
t69 = t75 * t87 + t98;
t68 = -g(3) * t87 + t91 * t84;
t65 = t70 * t86 + t71 * t83 + pkin(2);
t64 = t65 * t87 + t98;
t63 = -g(1) * (t85 * t76 + t87 * t94) - g(2) * (t77 * t96 - t95) - t77 * t99;
t1 = [0, -t91, g(1) * t85 - g(2) * t88, 0, 0, 0, 0, 0, -t91 * t87 - t99, t68, 0, 0, 0, 0, 0, -g(1) * (t85 * t83 + t87 * t92) - g(2) * (t86 * t96 - t93) - t86 * t99, -g(1) * (t85 * t86 - t87 * t93) - g(2) * (-t83 * t96 - t92) + t83 * t99, t63, t89, -t68, -g(1) * (t69 * t88 + t74 * t85) - g(2) * (t69 * t85 - t74 * t88) - g(3) * (t84 * t75 + t97), t63, -t68, -t89, -g(1) * (t64 * t88 + t90 * t85) - g(2) * (t64 * t85 - t90 * t88) - g(3) * (t65 * t84 + t97);];
U_reg = t1;
