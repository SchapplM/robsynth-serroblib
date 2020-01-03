% Calculate inertial parameters regressor of potential energy for
% S5RPRRP11
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:31
% EndTime: 2019-12-31 18:54:31
% DurationCPUTime: 0.11s
% Computational Cost: add. (129->49), mult. (156->61), div. (0->0), fcn. (155->8), ass. (0->34)
t71 = pkin(8) + qJ(3);
t67 = sin(t71);
t68 = cos(t71);
t95 = pkin(3) * t68 + pkin(7) * t67;
t94 = g(3) * pkin(5);
t91 = g(3) * t67;
t72 = sin(pkin(8));
t90 = t72 * pkin(2) + pkin(5);
t75 = sin(qJ(4));
t76 = sin(qJ(1));
t89 = t76 * t75;
t77 = cos(qJ(4));
t88 = t76 * t77;
t78 = cos(qJ(1));
t87 = t78 * t75;
t86 = t78 * t77;
t73 = cos(pkin(8));
t64 = t73 * pkin(2) + pkin(1);
t74 = -pkin(6) - qJ(2);
t85 = t76 * t64 + t78 * t74;
t84 = t78 * t64 - t76 * t74;
t83 = t95 * t76 + t85;
t82 = g(1) * t78 + g(2) * t76;
t81 = t67 * pkin(3) - t68 * pkin(7) + t90;
t80 = t95 * t78 + t84;
t51 = t68 * t89 + t86;
t53 = t68 * t87 - t88;
t79 = g(1) * t53 + g(2) * t51 + t75 * t91;
t57 = g(1) * t76 - g(2) * t78;
t54 = t68 * t86 + t89;
t52 = t68 * t88 - t87;
t50 = -g(3) * t68 + t82 * t67;
t49 = -g(1) * t54 - g(2) * t52 - t77 * t91;
t1 = [0, 0, 0, 0, 0, 0, -t82, t57, -g(3), -t94, 0, 0, 0, 0, 0, 0, -g(3) * t72 - t82 * t73, -g(3) * t73 + t82 * t72, -t57, -g(1) * (t78 * pkin(1) + t76 * qJ(2)) - g(2) * (t76 * pkin(1) - t78 * qJ(2)) - t94, 0, 0, 0, 0, 0, 0, -t82 * t68 - t91, t50, -t57, -g(1) * t84 - g(2) * t85 - g(3) * t90, 0, 0, 0, 0, 0, 0, t49, t79, -t50, -g(1) * t80 - g(2) * t83 - g(3) * t81, 0, 0, 0, 0, 0, 0, t49, -t50, -t79, -g(1) * (t54 * pkin(4) + t53 * qJ(5) + t80) - g(2) * (t52 * pkin(4) + t51 * qJ(5) + t83) - g(3) * ((pkin(4) * t77 + qJ(5) * t75) * t67 + t81);];
U_reg = t1;
