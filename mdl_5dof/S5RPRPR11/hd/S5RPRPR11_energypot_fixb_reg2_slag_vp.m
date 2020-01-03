% Calculate inertial parameters regressor of potential energy for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:54
% EndTime: 2019-12-31 18:27:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (117->46), mult. (133->56), div. (0->0), fcn. (126->8), ass. (0->29)
t77 = cos(qJ(1));
t70 = pkin(8) + qJ(3);
t66 = sin(t70);
t84 = qJ(4) * t66;
t67 = cos(t70);
t86 = t67 * t77;
t92 = pkin(3) * t86 + t77 * t84;
t75 = sin(qJ(1));
t81 = g(1) * t77 + g(2) * t75;
t91 = g(3) * pkin(5);
t71 = sin(pkin(8));
t88 = t71 * pkin(2) + pkin(5);
t87 = t67 * t75;
t72 = cos(pkin(8));
t64 = t72 * pkin(2) + pkin(1);
t73 = -pkin(6) - qJ(2);
t85 = t75 * t64 + t77 * t73;
t59 = t77 * t64;
t83 = -t75 * t73 + t59;
t82 = pkin(3) * t87 + t75 * t84 + t85;
t74 = sin(qJ(5));
t76 = cos(qJ(5));
t80 = t66 * t76 - t67 * t74;
t79 = t66 * t74 + t67 * t76;
t78 = t66 * pkin(3) - t67 * qJ(4) + t88;
t60 = g(1) * t75 - g(2) * t77;
t55 = -g(3) * t66 - t81 * t67;
t54 = -g(3) * t67 + t81 * t66;
t1 = [0, 0, 0, 0, 0, 0, -t81, t60, -g(3), -t91, 0, 0, 0, 0, 0, 0, -g(3) * t71 - t81 * t72, -g(3) * t72 + t81 * t71, -t60, -g(1) * (t77 * pkin(1) + t75 * qJ(2)) - g(2) * (t75 * pkin(1) - t77 * qJ(2)) - t91, 0, 0, 0, 0, 0, 0, t55, t54, -t60, -g(1) * t83 - g(2) * t85 - g(3) * t88, 0, 0, 0, 0, 0, 0, t55, -t60, -t54, -g(1) * (t83 + t92) - g(2) * t82 - g(3) * t78, 0, 0, 0, 0, 0, 0, -g(3) * t80 - t81 * t79, g(3) * t79 - t81 * t80, t60, -g(1) * (pkin(4) * t86 + t59 + (-pkin(7) - t73) * t75 + t92) - g(2) * (pkin(4) * t87 + t77 * pkin(7) + t82) - g(3) * (t66 * pkin(4) + t78);];
U_reg = t1;
