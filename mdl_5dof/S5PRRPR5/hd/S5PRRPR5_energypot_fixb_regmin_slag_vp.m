% Calculate minimal parameter regressor of potential energy for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:41
% EndTime: 2021-01-15 16:04:41
% DurationCPUTime: 0.12s
% Computational Cost: add. (113->56), mult. (206->97), div. (0->0), fcn. (252->12), ass. (0->33)
t78 = sin(pkin(9));
t79 = sin(pkin(5));
t99 = t78 * t79;
t80 = cos(pkin(9));
t98 = t79 * t80;
t84 = sin(qJ(3));
t97 = t79 * t84;
t85 = sin(qJ(2));
t96 = t79 * t85;
t87 = cos(qJ(3));
t95 = t79 * t87;
t88 = cos(qJ(2));
t94 = t79 * t88;
t81 = cos(pkin(5));
t93 = t81 * t84;
t92 = t81 * t85;
t91 = t81 * t88;
t69 = t78 * t85 - t80 * t91;
t71 = t78 * t91 + t80 * t85;
t89 = -g(1) * t71 - g(2) * t69 + g(3) * t94;
t86 = cos(qJ(5));
t83 = sin(qJ(5));
t82 = qJ(4) + pkin(7);
t77 = qJ(3) + pkin(10);
t76 = cos(t77);
t75 = sin(t77);
t74 = t87 * pkin(3) + pkin(2);
t70 = t78 * t88 + t80 * t92;
t68 = t78 * t92 - t80 * t88;
t67 = t81 * t75 + t76 * t96;
t66 = -t68 * t76 + t75 * t99;
t65 = t70 * t76 - t75 * t98;
t1 = [-g(3) * qJ(1), 0, g(1) * t68 - g(2) * t70 - g(3) * t96, -t89, 0, 0, 0, 0, 0, -g(1) * (-t68 * t87 + t78 * t97) - g(2) * (t70 * t87 - t80 * t97) - g(3) * (t85 * t95 + t93), -g(1) * (t68 * t84 + t78 * t95) - g(2) * (-t70 * t84 - t80 * t95) - g(3) * (t81 * t87 - t84 * t96), -g(1) * t66 - g(2) * t65 - g(3) * t67, -g(1) * (t68 * t75 + t76 * t99) - g(2) * (-t70 * t75 - t76 * t98) - g(3) * (-t75 * t96 + t81 * t76), t89, -g(1) * (t80 * pkin(1) - t68 * t74 + t71 * t82) - g(2) * (t78 * pkin(1) + t69 * t82 + t70 * t74) - g(3) * (pkin(3) * t93 + t81 * pkin(6) + qJ(1)) + (-g(3) * (t74 * t85 - t82 * t88) + (-g(1) * t78 + g(2) * t80) * (pkin(3) * t84 + pkin(6))) * t79, 0, 0, 0, 0, 0, -g(1) * (t66 * t86 + t71 * t83) - g(2) * (t65 * t86 + t69 * t83) - g(3) * (t67 * t86 - t83 * t94), -g(1) * (-t66 * t83 + t71 * t86) - g(2) * (-t65 * t83 + t69 * t86) - g(3) * (-t67 * t83 - t86 * t94);];
U_reg = t1;
