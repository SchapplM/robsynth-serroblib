% Calculate inertial parameters regressor of potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:06:51
% EndTime: 2019-12-05 18:06:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (122->64), mult. (169->83), div. (0->0), fcn. (168->8), ass. (0->38)
t97 = g(1) * pkin(5);
t79 = -pkin(7) - pkin(6);
t73 = sin(pkin(8));
t96 = g(1) * t73;
t76 = sin(qJ(1));
t95 = g(2) * t76;
t75 = sin(qJ(3));
t94 = t75 * pkin(3);
t77 = cos(qJ(3));
t64 = t77 * pkin(3) + pkin(2);
t71 = -qJ(5) + t79;
t93 = t71 * t73;
t92 = t73 * t79;
t74 = cos(pkin(8));
t78 = cos(qJ(1));
t91 = t74 * t78;
t72 = qJ(3) + qJ(4);
t65 = sin(t72);
t90 = t76 * t65;
t66 = cos(t72);
t89 = t76 * t66;
t88 = t76 * t75;
t87 = t76 * t77;
t86 = t78 * t65;
t85 = t78 * t66;
t84 = t78 * t75;
t83 = t78 * t77;
t82 = t78 * pkin(1) + t76 * qJ(2);
t81 = pkin(2) * t74 + pkin(6) * t73;
t80 = -g(3) * t78 + t95;
t68 = t78 * qJ(2);
t63 = g(2) * t78 + g(3) * t76;
t62 = pkin(4) * t65 + t94;
t61 = pkin(4) * t66 + t64;
t60 = g(1) * t74 + t80 * t73;
t59 = -t66 * t96 - g(2) * (-t74 * t89 + t86) - g(3) * (t74 * t85 + t90);
t58 = t65 * t96 - g(2) * (t74 * t90 + t85) - g(3) * (-t74 * t86 + t89);
t1 = [0, 0, 0, 0, 0, 0, t80, t63, -g(1), -t97, 0, 0, 0, 0, 0, 0, t80 * t74 - t96, -t60, -t63, -t97 - g(2) * (-t76 * pkin(1) + t68) - g(3) * t82, 0, 0, 0, 0, 0, 0, -t77 * t96 - g(2) * (-t74 * t87 + t84) - g(3) * (t74 * t83 + t88), t75 * t96 - g(2) * (t74 * t88 + t83) - g(3) * (-t74 * t84 + t87), t60, -g(1) * (t73 * pkin(2) - t74 * pkin(6) + pkin(5)) - g(2) * t68 - g(3) * (t81 * t78 + t82) - (-pkin(1) - t81) * t95, 0, 0, 0, 0, 0, 0, t59, t58, t60, -g(1) * (t73 * t64 + t74 * t79 + pkin(5)) - g(2) * (pkin(3) * t84 + t68) - g(3) * (t64 * t91 - t78 * t92 + t82) + (-g(2) * (-t64 * t74 - pkin(1) + t92) - g(3) * t94) * t76, 0, 0, 0, 0, 0, 0, t59, t58, t60, -g(1) * (t73 * t61 + t74 * t71 + pkin(5)) - g(2) * (t78 * t62 + t68) - g(3) * (t61 * t91 - t78 * t93 + t82) + (-g(2) * (-t61 * t74 - pkin(1) + t93) - g(3) * t62) * t76;];
U_reg = t1;
