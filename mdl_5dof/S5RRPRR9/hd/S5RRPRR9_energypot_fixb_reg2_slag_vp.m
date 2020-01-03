% Calculate inertial parameters regressor of potential energy for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:57
% EndTime: 2019-12-31 20:21:57
% DurationCPUTime: 0.18s
% Computational Cost: add. (131->61), mult. (143->79), div. (0->0), fcn. (138->10), ass. (0->36)
t77 = cos(qJ(4));
t63 = t77 * pkin(4) + pkin(3);
t71 = qJ(2) + pkin(9);
t65 = sin(t71);
t66 = cos(t71);
t80 = -pkin(8) - pkin(7);
t98 = t63 * t66 - t65 * t80;
t97 = g(3) * pkin(5);
t96 = g(3) * t65;
t75 = sin(qJ(2));
t95 = t75 * pkin(2) + pkin(5);
t72 = qJ(4) + qJ(5);
t67 = sin(t72);
t76 = sin(qJ(1));
t92 = t76 * t67;
t68 = cos(t72);
t91 = t76 * t68;
t74 = sin(qJ(4));
t90 = t76 * t74;
t89 = t76 * t77;
t79 = cos(qJ(1));
t88 = t79 * t67;
t87 = t79 * t68;
t86 = t79 * t74;
t85 = t79 * t77;
t78 = cos(qJ(2));
t64 = t78 * pkin(2) + pkin(1);
t73 = -pkin(6) - qJ(3);
t84 = t76 * t64 + t79 * t73;
t60 = t79 * t64;
t83 = -t76 * t73 + t60;
t82 = pkin(3) * t66 + pkin(7) * t65;
t81 = g(1) * t79 + g(2) * t76;
t58 = g(1) * t76 - g(2) * t79;
t57 = -g(3) * t66 + t65 * t81;
t1 = [0, 0, 0, 0, 0, 0, -t81, t58, -g(3), -t97, 0, 0, 0, 0, 0, 0, -g(3) * t75 - t78 * t81, -g(3) * t78 + t75 * t81, -t58, -g(1) * (t79 * pkin(1) + t76 * pkin(6)) - g(2) * (t76 * pkin(1) - t79 * pkin(6)) - t97, 0, 0, 0, 0, 0, 0, -t66 * t81 - t96, t57, -t58, -g(1) * t83 - g(2) * t84 - g(3) * t95, 0, 0, 0, 0, 0, 0, -g(1) * (t66 * t85 + t90) - g(2) * (t66 * t89 - t86) - t77 * t96, -g(1) * (-t66 * t86 + t89) - g(2) * (-t66 * t90 - t85) + t74 * t96, -t57, -g(1) * (t79 * t82 + t83) - g(2) * (t76 * t82 + t84) - g(3) * (t65 * pkin(3) - t66 * pkin(7) + t95), 0, 0, 0, 0, 0, 0, -g(1) * (t66 * t87 + t92) - g(2) * (t66 * t91 - t88) - t68 * t96, -g(1) * (-t66 * t88 + t91) - g(2) * (-t66 * t92 - t87) + t67 * t96, -t57, -g(1) * (t79 * t98 + t60) - g(2) * (-pkin(4) * t86 + t84) - g(3) * (t65 * t63 + t66 * t80 + t95) + (-g(1) * (pkin(4) * t74 - t73) - g(2) * t98) * t76;];
U_reg = t1;
