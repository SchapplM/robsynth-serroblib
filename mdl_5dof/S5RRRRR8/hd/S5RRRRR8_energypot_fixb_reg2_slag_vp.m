% Calculate inertial parameters regressor of potential energy for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:46
% EndTime: 2019-12-31 22:25:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (131->61), mult. (143->79), div. (0->0), fcn. (138->10), ass. (0->36)
t76 = cos(qJ(4));
t62 = t76 * pkin(4) + pkin(3);
t72 = qJ(2) + qJ(3);
t66 = sin(t72);
t68 = cos(t72);
t79 = -pkin(9) - pkin(8);
t98 = t62 * t68 - t66 * t79;
t97 = g(3) * pkin(5);
t96 = g(3) * t66;
t74 = sin(qJ(2));
t95 = t74 * pkin(2) + pkin(5);
t71 = qJ(4) + qJ(5);
t65 = sin(t71);
t75 = sin(qJ(1));
t92 = t75 * t65;
t67 = cos(t71);
t91 = t75 * t67;
t73 = sin(qJ(4));
t90 = t75 * t73;
t89 = t75 * t76;
t78 = cos(qJ(1));
t88 = t78 * t65;
t87 = t78 * t67;
t86 = t78 * t73;
t85 = t78 * t76;
t77 = cos(qJ(2));
t63 = t77 * pkin(2) + pkin(1);
t80 = -pkin(7) - pkin(6);
t84 = t75 * t63 + t78 * t80;
t60 = t78 * t63;
t83 = -t75 * t80 + t60;
t82 = pkin(3) * t68 + pkin(8) * t66;
t81 = g(1) * t78 + g(2) * t75;
t58 = g(1) * t75 - g(2) * t78;
t57 = -g(3) * t68 + t81 * t66;
t1 = [0, 0, 0, 0, 0, 0, -t81, t58, -g(3), -t97, 0, 0, 0, 0, 0, 0, -g(3) * t74 - t81 * t77, -g(3) * t77 + t81 * t74, -t58, -g(1) * (t78 * pkin(1) + t75 * pkin(6)) - g(2) * (t75 * pkin(1) - t78 * pkin(6)) - t97, 0, 0, 0, 0, 0, 0, -t81 * t68 - t96, t57, -t58, -g(1) * t83 - g(2) * t84 - g(3) * t95, 0, 0, 0, 0, 0, 0, -g(1) * (t68 * t85 + t90) - g(2) * (t68 * t89 - t86) - t76 * t96, -g(1) * (-t68 * t86 + t89) - g(2) * (-t68 * t90 - t85) + t73 * t96, -t57, -g(1) * (t82 * t78 + t83) - g(2) * (t82 * t75 + t84) - g(3) * (t66 * pkin(3) - t68 * pkin(8) + t95), 0, 0, 0, 0, 0, 0, -g(1) * (t68 * t87 + t92) - g(2) * (t68 * t91 - t88) - t67 * t96, -g(1) * (-t68 * t88 + t91) - g(2) * (-t68 * t92 - t87) + t65 * t96, -t57, -g(1) * (t78 * t98 + t60) - g(2) * (-pkin(4) * t86 + t84) - g(3) * (t66 * t62 + t68 * t79 + t95) + (-g(1) * (pkin(4) * t73 - t80) - g(2) * t98) * t75;];
U_reg = t1;
