% Calculate inertial parameters regressor of potential energy for
% S5RRRPP8
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:23
% EndTime: 2019-12-31 21:09:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (101->53), mult. (205->62), div. (0->0), fcn. (214->6), ass. (0->33)
t71 = sin(qJ(2));
t73 = cos(qJ(3));
t90 = t71 * t73;
t70 = sin(qJ(3));
t92 = t70 * t71;
t95 = pkin(3) * t90 + qJ(4) * t92;
t94 = g(3) * pkin(5);
t93 = t71 * pkin(2) + pkin(5);
t72 = sin(qJ(1));
t91 = t71 * t72;
t75 = cos(qJ(1));
t89 = t71 * t75;
t74 = cos(qJ(2));
t88 = t72 * t74;
t87 = t75 * t70;
t86 = t75 * t73;
t85 = t75 * pkin(1) + t72 * pkin(6);
t84 = t72 * pkin(1) - t75 * pkin(6);
t83 = t75 * t74 * pkin(2) + pkin(7) * t89 + t85;
t82 = -t74 * pkin(7) + t93;
t81 = g(1) * t75 + g(2) * t72;
t80 = pkin(2) * t88 + pkin(7) * t91 + t84;
t54 = -t72 * t73 + t74 * t87;
t55 = t72 * t70 + t74 * t86;
t79 = t55 * pkin(3) + t54 * qJ(4) + t83;
t52 = t70 * t88 + t86;
t78 = g(1) * t54 + g(2) * t52 + g(3) * t92;
t53 = t73 * t88 - t87;
t77 = g(1) * t55 + g(2) * t53 + g(3) * t90;
t76 = t53 * pkin(3) + t52 * qJ(4) + t80;
t56 = g(1) * t72 - g(2) * t75;
t49 = -g(3) * t74 + t81 * t71;
t1 = [0, 0, 0, 0, 0, 0, -t81, t56, -g(3), -t94, 0, 0, 0, 0, 0, 0, -g(3) * t71 - t81 * t74, t49, -t56, -g(1) * t85 - g(2) * t84 - t94, 0, 0, 0, 0, 0, 0, -t77, t78, -t49, -g(1) * t83 - g(2) * t80 - g(3) * t82, 0, 0, 0, 0, 0, 0, -t49, t77, -t78, -g(1) * t79 - g(2) * t76 - g(3) * (t82 + t95), 0, 0, 0, 0, 0, 0, -t49, -t78, -t77, -g(1) * (pkin(4) * t89 + t55 * qJ(5) + t79) - g(2) * (pkin(4) * t91 + t53 * qJ(5) + t76) - g(3) * (qJ(5) * t90 + (-pkin(4) - pkin(7)) * t74 + t93 + t95);];
U_reg = t1;
