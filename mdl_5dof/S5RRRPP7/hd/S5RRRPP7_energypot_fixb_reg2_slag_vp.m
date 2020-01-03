% Calculate inertial parameters regressor of potential energy for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:53
% EndTime: 2019-12-31 21:05:53
% DurationCPUTime: 0.12s
% Computational Cost: add. (101->50), mult. (205->61), div. (0->0), fcn. (214->6), ass. (0->33)
t70 = sin(qJ(2));
t72 = cos(qJ(3));
t88 = t70 * t72;
t69 = sin(qJ(3));
t89 = t69 * t70;
t93 = pkin(3) * t88 + qJ(4) * t89;
t92 = g(3) * pkin(5);
t91 = pkin(7) * t70;
t90 = t70 * pkin(2) + pkin(5);
t71 = sin(qJ(1));
t73 = cos(qJ(2));
t87 = t71 * t73;
t74 = cos(qJ(1));
t86 = t74 * t69;
t85 = t74 * t72;
t84 = t74 * pkin(1) + t71 * pkin(6);
t83 = qJ(5) * t70;
t82 = t71 * pkin(1) - t74 * pkin(6);
t81 = t84 + (pkin(2) * t73 + t91) * t74;
t80 = -t73 * pkin(7) + t90;
t79 = g(1) * t74 + g(2) * t71;
t78 = pkin(2) * t87 + t71 * t91 + t82;
t54 = -t71 * t72 + t73 * t86;
t55 = t71 * t69 + t73 * t85;
t77 = t55 * pkin(3) + t54 * qJ(4) + t81;
t52 = t69 * t87 + t85;
t76 = g(1) * t54 + g(2) * t52 + g(3) * t89;
t53 = t72 * t87 - t86;
t75 = t53 * pkin(3) + t52 * qJ(4) + t78;
t56 = g(1) * t71 - g(2) * t74;
t49 = -g(3) * t73 + t79 * t70;
t48 = -g(1) * t55 - g(2) * t53 - g(3) * t88;
t1 = [0, 0, 0, 0, 0, 0, -t79, t56, -g(3), -t92, 0, 0, 0, 0, 0, 0, -g(3) * t70 - t79 * t73, t49, -t56, -g(1) * t84 - g(2) * t82 - t92, 0, 0, 0, 0, 0, 0, t48, t76, -t49, -g(1) * t81 - g(2) * t78 - g(3) * t80, 0, 0, 0, 0, 0, 0, t48, -t49, -t76, -g(1) * t77 - g(2) * t75 - g(3) * (t80 + t93), 0, 0, 0, 0, 0, 0, t48, -t76, t49, -g(1) * (t55 * pkin(4) - t74 * t83 + t77) - g(2) * (t53 * pkin(4) - t71 * t83 + t75) - g(3) * (pkin(4) * t88 + (-pkin(7) + qJ(5)) * t73 + t90 + t93);];
U_reg = t1;
