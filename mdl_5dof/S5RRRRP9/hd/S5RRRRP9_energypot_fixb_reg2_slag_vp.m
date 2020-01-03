% Calculate inertial parameters regressor of potential energy for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:03
% EndTime: 2019-12-31 22:06:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (129->61), mult. (184->77), div. (0->0), fcn. (187->8), ass. (0->37)
t99 = g(3) * pkin(5);
t77 = sin(qJ(2));
t98 = g(3) * t77;
t82 = -pkin(8) - pkin(7);
t97 = t77 * t82;
t76 = sin(qJ(3));
t78 = sin(qJ(1));
t96 = t78 * t76;
t80 = cos(qJ(2));
t95 = t78 * t80;
t75 = qJ(3) + qJ(4);
t69 = sin(t75);
t81 = cos(qJ(1));
t94 = t81 * t69;
t70 = cos(t75);
t93 = t81 * t70;
t92 = t81 * t76;
t79 = cos(qJ(3));
t91 = t81 * t79;
t90 = t81 * pkin(1) + t78 * pkin(6);
t67 = t79 * pkin(3) + pkin(2);
t89 = t77 * t67 + t80 * t82 + pkin(5);
t72 = t78 * pkin(1);
t88 = -t81 * pkin(6) + t72;
t87 = pkin(2) * t80 + pkin(7) * t77;
t86 = g(1) * t81 + g(2) * t78;
t85 = pkin(3) * t96 + t90 + (t67 * t80 - t97) * t81;
t56 = t69 * t95 + t93;
t58 = -t78 * t70 + t80 * t94;
t84 = g(1) * t58 + g(2) * t56 + t69 * t98;
t83 = -t78 * t97 + t67 * t95 + t72 + (-pkin(3) * t76 - pkin(6)) * t81;
t63 = g(1) * t78 - g(2) * t81;
t60 = -g(3) * t80 + t86 * t77;
t59 = t78 * t69 + t80 * t93;
t57 = t70 * t95 - t94;
t55 = -g(1) * t59 - g(2) * t57 - t70 * t98;
t1 = [0, 0, 0, 0, 0, 0, -t86, t63, -g(3), -t99, 0, 0, 0, 0, 0, 0, -t86 * t80 - t98, t60, -t63, -g(1) * t90 - g(2) * t88 - t99, 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t91 + t96) - g(2) * (t79 * t95 - t92) - t79 * t98, -g(1) * (t78 * t79 - t80 * t92) - g(2) * (-t76 * t95 - t91) + t76 * t98, -t60, -g(1) * (t87 * t81 + t90) - g(2) * (t87 * t78 + t88) - g(3) * (t77 * pkin(2) - t80 * pkin(7) + pkin(5)), 0, 0, 0, 0, 0, 0, t55, t84, -t60, -g(1) * t85 - g(2) * t83 - g(3) * t89, 0, 0, 0, 0, 0, 0, t55, -t60, -t84, -g(1) * (t59 * pkin(4) + t58 * qJ(5) + t85) - g(2) * (t57 * pkin(4) + t56 * qJ(5) + t83) - g(3) * ((pkin(4) * t70 + qJ(5) * t69) * t77 + t89);];
U_reg = t1;
