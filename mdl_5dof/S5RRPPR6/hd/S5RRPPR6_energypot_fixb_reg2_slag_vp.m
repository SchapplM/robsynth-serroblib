% Calculate inertial parameters regressor of potential energy for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:04
% EndTime: 2019-12-31 19:33:04
% DurationCPUTime: 0.18s
% Computational Cost: add. (131->61), mult. (143->79), div. (0->0), fcn. (138->10), ass. (0->36)
t74 = cos(pkin(9));
t62 = t74 * pkin(4) + pkin(3);
t72 = qJ(2) + pkin(8);
t66 = sin(t72);
t68 = cos(t72);
t75 = -pkin(7) - qJ(4);
t98 = t62 * t68 - t66 * t75;
t97 = g(3) * pkin(5);
t96 = g(3) * t66;
t77 = sin(qJ(2));
t95 = t77 * pkin(2) + pkin(5);
t71 = pkin(9) + qJ(5);
t65 = sin(t71);
t78 = sin(qJ(1));
t92 = t78 * t65;
t67 = cos(t71);
t91 = t78 * t67;
t73 = sin(pkin(9));
t90 = t78 * t73;
t89 = t78 * t74;
t80 = cos(qJ(1));
t88 = t80 * t65;
t87 = t80 * t67;
t86 = t80 * t73;
t85 = t80 * t74;
t79 = cos(qJ(2));
t64 = t79 * pkin(2) + pkin(1);
t76 = -pkin(6) - qJ(3);
t84 = t78 * t64 + t80 * t76;
t60 = t80 * t64;
t83 = -t78 * t76 + t60;
t82 = g(1) * t80 + g(2) * t78;
t81 = pkin(3) * t68 + qJ(4) * t66;
t58 = g(1) * t78 - g(2) * t80;
t57 = -g(3) * t68 + t82 * t66;
t1 = [0, 0, 0, 0, 0, 0, -t82, t58, -g(3), -t97, 0, 0, 0, 0, 0, 0, -g(3) * t77 - t82 * t79, -g(3) * t79 + t82 * t77, -t58, -g(1) * (t80 * pkin(1) + t78 * pkin(6)) - g(2) * (t78 * pkin(1) - t80 * pkin(6)) - t97, 0, 0, 0, 0, 0, 0, -t82 * t68 - t96, t57, -t58, -g(1) * t83 - g(2) * t84 - g(3) * t95, 0, 0, 0, 0, 0, 0, -g(1) * (t68 * t85 + t90) - g(2) * (t68 * t89 - t86) - t74 * t96, -g(1) * (-t68 * t86 + t89) - g(2) * (-t68 * t90 - t85) + t73 * t96, -t57, -g(1) * (t81 * t80 + t83) - g(2) * (t81 * t78 + t84) - g(3) * (t66 * pkin(3) - t68 * qJ(4) + t95), 0, 0, 0, 0, 0, 0, -g(1) * (t68 * t87 + t92) - g(2) * (t68 * t91 - t88) - t67 * t96, -g(1) * (-t68 * t88 + t91) - g(2) * (-t68 * t92 - t87) + t65 * t96, -t57, -g(1) * (t80 * t98 + t60) - g(2) * (-pkin(4) * t86 + t84) - g(3) * (t66 * t62 + t68 * t75 + t95) + (-g(1) * (pkin(4) * t73 - t76) - g(2) * t98) * t78;];
U_reg = t1;
