% Calculate inertial parameters regressor of potential energy for
% S5RRRRP8
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:07
% EndTime: 2019-12-31 22:02:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (122->66), mult. (169->83), div. (0->0), fcn. (168->8), ass. (0->35)
t93 = g(3) * pkin(5);
t78 = -pkin(8) - pkin(7);
t73 = sin(qJ(2));
t92 = g(3) * t73;
t72 = sin(qJ(3));
t91 = t72 * pkin(3);
t75 = cos(qJ(3));
t62 = t75 * pkin(3) + pkin(2);
t70 = -qJ(5) + t78;
t90 = t70 * t73;
t89 = t73 * t78;
t74 = sin(qJ(1));
t88 = t74 * t72;
t76 = cos(qJ(2));
t87 = t74 * t76;
t71 = qJ(3) + qJ(4);
t63 = sin(t71);
t77 = cos(qJ(1));
t86 = t77 * t63;
t64 = cos(t71);
t85 = t77 * t64;
t84 = t77 * t72;
t83 = t77 * t75;
t82 = t77 * pkin(1) + t74 * pkin(6);
t66 = t74 * pkin(1);
t81 = -t77 * pkin(6) + t66;
t80 = pkin(2) * t76 + pkin(7) * t73;
t79 = g(1) * t77 + g(2) * t74;
t61 = g(1) * t74 - g(2) * t77;
t60 = pkin(4) * t63 + t91;
t59 = pkin(4) * t64 + t62;
t58 = -g(3) * t76 + t79 * t73;
t57 = -g(1) * (t74 * t63 + t76 * t85) - g(2) * (t64 * t87 - t86) - t64 * t92;
t56 = -g(1) * (t74 * t64 - t76 * t86) - g(2) * (-t63 * t87 - t85) + t63 * t92;
t1 = [0, 0, 0, 0, 0, 0, -t79, t61, -g(3), -t93, 0, 0, 0, 0, 0, 0, -t79 * t76 - t92, t58, -t61, -g(1) * t82 - g(2) * t81 - t93, 0, 0, 0, 0, 0, 0, -g(1) * (t76 * t83 + t88) - g(2) * (t75 * t87 - t84) - t75 * t92, -g(1) * (t74 * t75 - t76 * t84) - g(2) * (-t72 * t87 - t83) + t72 * t92, -t58, -g(1) * (t80 * t77 + t82) - g(2) * (t80 * t74 + t81) - g(3) * (t73 * pkin(2) - t76 * pkin(7) + pkin(5)), 0, 0, 0, 0, 0, 0, t57, t56, -t58, -g(1) * (pkin(3) * t88 + t82) - g(2) * (t62 * t87 - t74 * t89 + t66) - g(3) * (t73 * t62 + t76 * t78 + pkin(5)) + (-g(1) * (t62 * t76 - t89) - g(2) * (-pkin(6) - t91)) * t77, 0, 0, 0, 0, 0, 0, t57, t56, -t58, -g(1) * (t74 * t60 + t82) - g(2) * (t59 * t87 - t74 * t90 + t66) - g(3) * (t73 * t59 + t76 * t70 + pkin(5)) + (-g(1) * (t59 * t76 - t90) - g(2) * (-pkin(6) - t60)) * t77;];
U_reg = t1;
