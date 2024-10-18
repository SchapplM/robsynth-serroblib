% Calculate inertial parameters regressor of potential energy for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR14_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:30
% EndTime: 2024-09-27 18:43:31
% DurationCPUTime: 0.05s
% Computational Cost: add. (223->80), mult. (167->110), div. (0->0), fcn. (154->22), ass. (0->49)
t108 = pkin(7) + pkin(6);
t100 = pkin(8) + pkin(9);
t93 = sin(pkin(5));
t107 = pkin(8) * t93;
t98 = cos(qJ(3));
t79 = t98 * pkin(3) + pkin(2);
t96 = sin(qJ(3));
t106 = t93 * t96;
t94 = cos(pkin(5));
t105 = t94 * t96;
t104 = t94 * t98;
t91 = qJ(3) + qJ(4);
t81 = pkin(5) - t91;
t80 = pkin(5) + t91;
t92 = qJ(1) + qJ(2);
t83 = sin(t92);
t85 = cos(t92);
t103 = g(1) * t83 - g(2) * t85;
t97 = sin(qJ(1));
t99 = cos(qJ(1));
t102 = -g(1) * t99 - g(2) * t97;
t101 = pkin(4) * sin(qJ(4)) * t98 + (cos(qJ(4)) * pkin(4) + pkin(3)) * t96;
t90 = pkin(10) + t100;
t89 = t99 * pkin(1);
t87 = t97 * pkin(1);
t86 = qJ(5) + t91;
t84 = cos(t91);
t82 = sin(t91);
t77 = cos(t86);
t76 = sin(t86);
t75 = -qJ(5) + t81;
t74 = qJ(5) + t80;
t73 = cos(t80);
t72 = sin(t81);
t71 = cos(t74);
t70 = sin(t75);
t69 = cos(t81) / 0.2e1;
t68 = sin(t80) / 0.2e1;
t67 = cos(t75) / 0.2e1;
t66 = sin(t74) / 0.2e1;
t65 = pkin(4) * t84 + t79;
t64 = pkin(3) * t105 - t100 * t93;
t63 = t69 + t73 / 0.2e1;
t62 = t68 - t72 / 0.2e1;
t61 = t67 + t71 / 0.2e1;
t60 = t66 - t70 / 0.2e1;
t59 = -g(3) * t94 - t103 * t93;
t58 = t101 * t94 - t93 * t90;
t1 = [0, 0, 0, 0, 0, 0, t102, g(1) * t97 - g(2) * t99, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -g(1) * t85 - g(2) * t83, t103, -g(3), pkin(1) * t102 - g(3) * t108, 0, 0, 0, 0, 0, 0, -g(1) * (-t105 * t83 + t85 * t98) - g(2) * (t105 * t85 + t83 * t98) - g(3) * t106, -g(1) * (-t104 * t83 - t85 * t96) - g(2) * (t104 * t85 - t83 * t96) - g(3) * t93 * t98, t59, -g(1) * (pkin(2) * t85 + t107 * t83 + t89) - g(2) * (pkin(2) * t83 - t107 * t85 + t87) - g(3) * (pkin(8) * t94 + t108), 0, 0, 0, 0, 0, 0, -g(1) * (-t62 * t83 + t84 * t85) - g(2) * (t62 * t85 + t83 * t84) - g(3) * (t69 - t73 / 0.2e1), -g(1) * (-t63 * t83 - t82 * t85) - g(2) * (t63 * t85 - t82 * t83) - g(3) * (t68 + t72 / 0.2e1), t59, -g(1) * (-t64 * t83 + t79 * t85 + t89) - g(2) * (t64 * t85 + t79 * t83 + t87) - g(3) * (pkin(3) * t106 + t100 * t94 + t108), 0, 0, 0, 0, 0, 0, -g(1) * (-t60 * t83 + t77 * t85) - g(2) * (t60 * t85 + t77 * t83) - g(3) * (t67 - t71 / 0.2e1), -g(1) * (-t61 * t83 - t76 * t85) - g(2) * (t61 * t85 - t76 * t83) - g(3) * (t66 + t70 / 0.2e1), t59, -g(1) * (-t58 * t83 + t65 * t85 + t89) - g(2) * (t58 * t85 + t65 * t83 + t87) - g(3) * (t101 * t93 + t94 * t90 + t108);];
U_reg = t1;
