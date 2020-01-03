% Calculate inertial parameters regressor of potential energy for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:21
% EndTime: 2019-12-31 19:36:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (113->52), mult. (130->62), div. (0->0), fcn. (121->8), ass. (0->32)
t78 = cos(qJ(1));
t71 = qJ(2) + pkin(8);
t67 = sin(t71);
t83 = qJ(4) * t67;
t68 = cos(t71);
t89 = t68 * t78;
t94 = pkin(3) * t89 + t78 * t83;
t93 = g(3) * pkin(5);
t92 = g(3) * t68;
t74 = sin(qJ(2));
t91 = t74 * pkin(2) + pkin(5);
t75 = sin(qJ(1));
t90 = t68 * t75;
t73 = sin(qJ(5));
t88 = t75 * t73;
t76 = cos(qJ(5));
t87 = t75 * t76;
t86 = t78 * t73;
t85 = t78 * t76;
t77 = cos(qJ(2));
t66 = t77 * pkin(2) + pkin(1);
t72 = -pkin(6) - qJ(3);
t84 = t75 * t66 + t78 * t72;
t63 = t78 * t66;
t82 = -t75 * t72 + t63;
t81 = pkin(3) * t90 + t75 * t83 + t84;
t80 = g(1) * t78 + g(2) * t75;
t79 = t67 * pkin(3) - t68 * qJ(4) + t91;
t59 = g(1) * t75 - g(2) * t78;
t56 = g(3) * t67 + t80 * t68;
t55 = t80 * t67 - t92;
t1 = [0, 0, 0, 0, 0, 0, -t80, t59, -g(3), -t93, 0, 0, 0, 0, 0, 0, -g(3) * t74 - t80 * t77, -g(3) * t77 + t80 * t74, -t59, -g(1) * (t78 * pkin(1) + t75 * pkin(6)) - g(2) * (t75 * pkin(1) - t78 * pkin(6)) - t93, 0, 0, 0, 0, 0, 0, -t56, t55, -t59, -g(1) * t82 - g(2) * t84 - g(3) * t91, 0, 0, 0, 0, 0, 0, -t59, t56, -t55, -g(1) * (t82 + t94) - g(2) * t81 - g(3) * t79, 0, 0, 0, 0, 0, 0, -g(1) * (t67 * t86 + t87) - g(2) * (t67 * t88 - t85) + t73 * t92, -g(1) * (t67 * t85 - t88) - g(2) * (t67 * t87 + t86) + t76 * t92, -t56, -g(1) * (pkin(7) * t89 + t63 + (pkin(4) - t72) * t75 + t94) - g(2) * (-t78 * pkin(4) + pkin(7) * t90 + t81) - g(3) * (t67 * pkin(7) + t79);];
U_reg = t1;
