% Calculate inertial parameters regressor of potential energy for
% S5RPRRR7
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:08
% EndTime: 2019-12-31 19:04:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (138->58), mult. (130->78), div. (0->0), fcn. (125->10), ass. (0->32)
t68 = cos(qJ(4));
t54 = t68 * pkin(4) + pkin(3);
t66 = sin(qJ(3));
t69 = cos(qJ(3));
t71 = -pkin(8) - pkin(7);
t87 = t54 * t69 - t66 * t71;
t64 = qJ(2) + pkin(5);
t86 = g(3) * t64;
t85 = g(3) * t66;
t62 = qJ(1) + pkin(9);
t55 = sin(t62);
t65 = sin(qJ(4));
t83 = t55 * t65;
t63 = qJ(4) + qJ(5);
t57 = sin(t63);
t82 = t57 * t69;
t58 = cos(t63);
t81 = t58 * t69;
t80 = t65 * t69;
t78 = t68 * t69;
t67 = sin(qJ(1));
t77 = t67 * pkin(1) + t55 * pkin(2);
t56 = cos(t62);
t70 = cos(qJ(1));
t76 = t70 * pkin(1) + t56 * pkin(2) + t55 * pkin(6);
t75 = -t56 * pkin(6) + t77;
t74 = pkin(3) * t69 + pkin(7) * t66;
t73 = g(1) * t56 + g(2) * t55;
t72 = -g(1) * t70 - g(2) * t67;
t50 = g(1) * t55 - g(2) * t56;
t49 = -g(3) * t69 + t73 * t66;
t1 = [0, 0, 0, 0, 0, 0, t72, g(1) * t67 - g(2) * t70, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t73, t50, -g(3), t72 * pkin(1) - t86, 0, 0, 0, 0, 0, 0, -t73 * t69 - t85, t49, -t50, -g(1) * t76 - g(2) * t75 - t86, 0, 0, 0, 0, 0, 0, -g(1) * (t56 * t78 + t83) - g(2) * (t55 * t78 - t56 * t65) - t68 * t85, -g(1) * (t55 * t68 - t56 * t80) - g(2) * (-t55 * t80 - t56 * t68) + t65 * t85, -t49, -g(1) * (t74 * t56 + t76) - g(2) * (t74 * t55 + t75) - g(3) * (t66 * pkin(3) - t69 * pkin(7) + t64), 0, 0, 0, 0, 0, 0, -g(1) * (t55 * t57 + t56 * t81) - g(2) * (t55 * t81 - t56 * t57) - t58 * t85, -g(1) * (t55 * t58 - t56 * t82) - g(2) * (-t55 * t82 - t56 * t58) + t57 * t85, -t49, -g(1) * (pkin(4) * t83 + t76) - g(2) * (t87 * t55 + t77) - g(3) * (t66 * t54 + t69 * t71 + t64) + (-g(1) * t87 - g(2) * (-pkin(4) * t65 - pkin(6))) * t56;];
U_reg = t1;
