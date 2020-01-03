% Calculate inertial parameters regressor of potential energy for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:58
% EndTime: 2019-12-31 19:07:58
% DurationCPUTime: 0.10s
% Computational Cost: add. (130->50), mult. (119->64), div. (0->0), fcn. (110->10), ass. (0->31)
t88 = g(3) * pkin(5);
t69 = pkin(9) + qJ(3);
t64 = qJ(4) + t69;
t59 = sin(t64);
t87 = g(3) * t59;
t70 = sin(pkin(9));
t86 = t70 * pkin(2) + pkin(5);
t71 = cos(pkin(9));
t61 = t71 * pkin(2) + pkin(1);
t73 = sin(qJ(5));
t74 = sin(qJ(1));
t85 = t74 * t73;
t75 = cos(qJ(5));
t84 = t74 * t75;
t76 = cos(qJ(1));
t83 = t76 * t73;
t82 = t76 * t75;
t72 = -pkin(6) - qJ(2);
t63 = cos(t69);
t55 = pkin(3) * t63 + t61;
t68 = -pkin(7) + t72;
t81 = t74 * t55 + t76 * t68;
t62 = sin(t69);
t80 = pkin(3) * t62 + t86;
t79 = t76 * t55 - t74 * t68;
t60 = cos(t64);
t78 = pkin(4) * t60 + pkin(8) * t59;
t77 = g(1) * t76 + g(2) * t74;
t56 = g(1) * t74 - g(2) * t76;
t52 = -g(3) * t60 + t77 * t59;
t1 = [0, 0, 0, 0, 0, 0, -t77, t56, -g(3), -t88, 0, 0, 0, 0, 0, 0, -g(3) * t70 - t77 * t71, -g(3) * t71 + t77 * t70, -t56, -g(1) * (t76 * pkin(1) + t74 * qJ(2)) - g(2) * (t74 * pkin(1) - t76 * qJ(2)) - t88, 0, 0, 0, 0, 0, 0, -g(3) * t62 - t77 * t63, -g(3) * t63 + t77 * t62, -t56, -g(1) * (t76 * t61 - t74 * t72) - g(2) * (t74 * t61 + t76 * t72) - g(3) * t86, 0, 0, 0, 0, 0, 0, -t77 * t60 - t87, t52, -t56, -g(1) * t79 - g(2) * t81 - g(3) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (t60 * t82 + t85) - g(2) * (t60 * t84 - t83) - t75 * t87, -g(1) * (-t60 * t83 + t84) - g(2) * (-t60 * t85 - t82) + t73 * t87, -t52, -g(1) * (t78 * t76 + t79) - g(2) * (t78 * t74 + t81) - g(3) * (t59 * pkin(4) - t60 * pkin(8) + t80);];
U_reg = t1;
