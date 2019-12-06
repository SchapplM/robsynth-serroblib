% Calculate inertial parameters regressor of potential energy for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:34
% EndTime: 2019-12-05 17:00:35
% DurationCPUTime: 0.16s
% Computational Cost: add. (202->66), mult. (470->99), div. (0->0), fcn. (575->10), ass. (0->48)
t83 = sin(pkin(5));
t111 = pkin(6) * t83;
t87 = sin(qJ(3));
t110 = t83 * t87;
t88 = sin(qJ(2));
t109 = t83 * t88;
t90 = cos(qJ(3));
t108 = t83 * t90;
t91 = cos(qJ(2));
t107 = t83 * t91;
t85 = cos(pkin(5));
t106 = t85 * t88;
t105 = t85 * t91;
t82 = sin(pkin(9));
t84 = cos(pkin(9));
t104 = t84 * pkin(1) + t82 * t111;
t103 = t85 * pkin(6) + qJ(1);
t102 = t82 * pkin(1) - t84 * t111;
t101 = g(1) * t82 - g(2) * t84;
t70 = t82 * t105 + t84 * t88;
t71 = -t82 * t106 + t84 * t91;
t100 = t71 * pkin(2) + t70 * pkin(7) + t104;
t99 = pkin(2) * t109 - pkin(7) * t107 + t103;
t69 = t84 * t106 + t82 * t91;
t58 = -t84 * t110 + t69 * t90;
t68 = -t84 * t105 + t82 * t88;
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t51 = t58 * t86 - t68 * t89;
t60 = t82 * t110 + t71 * t90;
t53 = t60 * t86 - t70 * t89;
t73 = t88 * t108 + t85 * t87;
t61 = t89 * t107 + t73 * t86;
t98 = g(1) * t53 + g(2) * t51 + g(3) * t61;
t57 = t84 * t108 + t69 * t87;
t59 = -t82 * t108 + t71 * t87;
t72 = t87 * t109 - t85 * t90;
t97 = g(1) * t59 + g(2) * t57 + g(3) * t72;
t96 = t69 * pkin(2) + t68 * pkin(7) + t102;
t95 = -g(1) * t70 - g(2) * t68 + g(3) * t107;
t94 = t60 * pkin(3) + t59 * pkin(8) + t100;
t93 = t73 * pkin(3) + t72 * pkin(8) + t99;
t92 = t58 * pkin(3) + t57 * pkin(8) + t96;
t62 = -t86 * t107 + t73 * t89;
t54 = t60 * t89 + t70 * t86;
t52 = t58 * t89 + t68 * t86;
t49 = -g(1) * t54 - g(2) * t52 - g(3) * t62;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t84 - g(2) * t82, t101, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t71 - g(2) * t69 - g(3) * t109, -t95, -g(3) * t85 - t101 * t83, -g(1) * t104 - g(2) * t102 - g(3) * t103, 0, 0, 0, 0, 0, 0, -g(1) * t60 - g(2) * t58 - g(3) * t73, t97, t95, -g(1) * t100 - g(2) * t96 - g(3) * t99, 0, 0, 0, 0, 0, 0, t49, t98, -t97, -g(1) * t94 - g(2) * t92 - g(3) * t93, 0, 0, 0, 0, 0, 0, t49, -t97, -t98, -g(1) * (t54 * pkin(4) + t53 * qJ(5) + t94) - g(2) * (t52 * pkin(4) + t51 * qJ(5) + t92) - g(3) * (t62 * pkin(4) + t61 * qJ(5) + t93);];
U_reg = t1;
