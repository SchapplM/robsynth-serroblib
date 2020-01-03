% Calculate inertial parameters regressor of potential energy for
% S5RPRRP10
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:03
% EndTime: 2019-12-31 18:52:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (121->52), mult. (143->65), div. (0->0), fcn. (138->8), ass. (0->31)
t71 = cos(qJ(4));
t59 = t71 * pkin(4) + pkin(3);
t64 = pkin(8) + qJ(3);
t60 = sin(t64);
t61 = cos(t64);
t67 = -qJ(5) - pkin(7);
t86 = t59 * t61 - t60 * t67;
t85 = g(3) * pkin(5);
t84 = g(3) * t60;
t65 = sin(pkin(8));
t83 = t65 * pkin(2) + pkin(5);
t69 = sin(qJ(4));
t70 = sin(qJ(1));
t80 = t70 * t69;
t79 = t70 * t71;
t72 = cos(qJ(1));
t78 = t72 * t69;
t77 = t72 * t71;
t66 = cos(pkin(8));
t57 = t66 * pkin(2) + pkin(1);
t68 = -pkin(6) - qJ(2);
t76 = t70 * t57 + t72 * t68;
t54 = t72 * t57;
t75 = -t70 * t68 + t54;
t74 = pkin(3) * t61 + pkin(7) * t60;
t73 = g(1) * t72 + g(2) * t70;
t55 = g(1) * t70 - g(2) * t72;
t52 = -g(3) * t61 + t73 * t60;
t51 = -g(1) * (t61 * t77 + t80) - g(2) * (t61 * t79 - t78) - t71 * t84;
t50 = -g(1) * (-t61 * t78 + t79) - g(2) * (-t61 * t80 - t77) + t69 * t84;
t1 = [0, 0, 0, 0, 0, 0, -t73, t55, -g(3), -t85, 0, 0, 0, 0, 0, 0, -g(3) * t65 - t73 * t66, -g(3) * t66 + t73 * t65, -t55, -g(1) * (t72 * pkin(1) + t70 * qJ(2)) - g(2) * (t70 * pkin(1) - t72 * qJ(2)) - t85, 0, 0, 0, 0, 0, 0, -t73 * t61 - t84, t52, -t55, -g(1) * t75 - g(2) * t76 - g(3) * t83, 0, 0, 0, 0, 0, 0, t51, t50, -t52, -g(1) * (t74 * t72 + t75) - g(2) * (t74 * t70 + t76) - g(3) * (t60 * pkin(3) - t61 * pkin(7) + t83), 0, 0, 0, 0, 0, 0, t51, t50, -t52, -g(1) * (t86 * t72 + t54) - g(2) * (-pkin(4) * t78 + t76) - g(3) * (t60 * t59 + t61 * t67 + t83) + (-g(1) * (pkin(4) * t69 - t68) - g(2) * t86) * t70;];
U_reg = t1;
