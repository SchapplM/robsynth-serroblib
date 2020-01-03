% Calculate inertial parameters regressor of potential energy for
% S5RPRPR12
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
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:14
% EndTime: 2019-12-31 18:30:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (131->61), mult. (143->79), div. (0->0), fcn. (138->10), ass. (0->36)
t71 = cos(pkin(9));
t58 = t71 * pkin(4) + pkin(3);
t68 = pkin(8) + qJ(3);
t62 = sin(t68);
t64 = cos(t68);
t73 = -pkin(7) - qJ(4);
t94 = t58 * t64 - t62 * t73;
t93 = g(3) * pkin(5);
t92 = g(3) * t62;
t70 = sin(pkin(8));
t91 = t70 * pkin(2) + pkin(5);
t67 = pkin(9) + qJ(5);
t61 = sin(t67);
t75 = sin(qJ(1));
t88 = t75 * t61;
t63 = cos(t67);
t87 = t75 * t63;
t69 = sin(pkin(9));
t86 = t75 * t69;
t85 = t75 * t71;
t76 = cos(qJ(1));
t84 = t76 * t61;
t83 = t76 * t63;
t82 = t76 * t69;
t81 = t76 * t71;
t72 = cos(pkin(8));
t59 = t72 * pkin(2) + pkin(1);
t74 = -pkin(6) - qJ(2);
t80 = t75 * t59 + t76 * t74;
t55 = t76 * t59;
t79 = -t75 * t74 + t55;
t78 = g(1) * t76 + g(2) * t75;
t77 = pkin(3) * t64 + qJ(4) * t62;
t56 = g(1) * t75 - g(2) * t76;
t53 = -g(3) * t64 + t78 * t62;
t1 = [0, 0, 0, 0, 0, 0, -t78, t56, -g(3), -t93, 0, 0, 0, 0, 0, 0, -g(3) * t70 - t78 * t72, -g(3) * t72 + t78 * t70, -t56, -g(1) * (t76 * pkin(1) + t75 * qJ(2)) - g(2) * (t75 * pkin(1) - t76 * qJ(2)) - t93, 0, 0, 0, 0, 0, 0, -t78 * t64 - t92, t53, -t56, -g(1) * t79 - g(2) * t80 - g(3) * t91, 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t81 + t86) - g(2) * (t64 * t85 - t82) - t71 * t92, -g(1) * (-t64 * t82 + t85) - g(2) * (-t64 * t86 - t81) + t69 * t92, -t53, -g(1) * (t77 * t76 + t79) - g(2) * (t77 * t75 + t80) - g(3) * (t62 * pkin(3) - t64 * qJ(4) + t91), 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t83 + t88) - g(2) * (t64 * t87 - t84) - t63 * t92, -g(1) * (-t64 * t84 + t87) - g(2) * (-t64 * t88 - t83) + t61 * t92, -t53, -g(1) * (t94 * t76 + t55) - g(2) * (-pkin(4) * t82 + t80) - g(3) * (t62 * t58 + t64 * t73 + t91) + (-g(1) * (pkin(4) * t69 - t74) - g(2) * t94) * t75;];
U_reg = t1;
