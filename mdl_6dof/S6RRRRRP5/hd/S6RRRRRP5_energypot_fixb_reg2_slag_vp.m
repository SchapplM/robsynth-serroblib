% Calculate inertial parameters regressor of potential energy for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:23:02
% EndTime: 2019-03-10 01:23:03
% DurationCPUTime: 0.18s
% Computational Cost: add. (211->90), mult. (228->114), div. (0->0), fcn. (228->10), ass. (0->40)
t116 = g(3) * pkin(6);
t103 = -pkin(9) - pkin(8);
t98 = sin(qJ(2));
t115 = g(3) * t98;
t97 = sin(qJ(3));
t89 = t97 * pkin(3);
t100 = cos(qJ(3));
t84 = t100 * pkin(3) + pkin(2);
t99 = sin(qJ(1));
t114 = t98 * t99;
t113 = t99 * t97;
t96 = qJ(3) + qJ(4);
t85 = sin(t96);
t78 = pkin(4) * t85 + t89;
t102 = cos(qJ(1));
t112 = t102 * pkin(1) + t99 * pkin(7);
t101 = cos(qJ(2));
t111 = t101 * t99;
t110 = t103 * t98;
t109 = t99 * t100;
t108 = t101 * t102;
t107 = t102 * t100;
t95 = -pkin(10) + t103;
t86 = cos(t96);
t77 = pkin(4) * t86 + t84;
t91 = t99 * pkin(1);
t106 = -t102 * pkin(7) + t91;
t105 = pkin(2) * t101 + pkin(8) * t98;
t104 = g(1) * t102 + g(2) * t99;
t88 = qJ(5) + t96;
t87 = -qJ(6) + t95;
t83 = cos(t88);
t82 = sin(t88);
t79 = g(1) * t99 - g(2) * t102;
t76 = pkin(5) * t82 + t78;
t75 = -g(3) * t101 + t104 * t98;
t74 = pkin(5) * t83 + t77;
t73 = -g(1) * (t83 * t108 + t99 * t82) - g(2) * (-t102 * t82 + t83 * t111) - t83 * t115;
t72 = -g(1) * (-t82 * t108 + t99 * t83) - g(2) * (-t102 * t83 - t82 * t111) + t82 * t115;
t1 = [0, 0, 0, 0, 0, 0, -t104, t79, -g(3), -t116, 0, 0, 0, 0, 0, 0, -t104 * t101 - t115, t75, -t79, -g(1) * t112 - g(2) * t106 - t116, 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t107 + t113) - g(2) * (t101 * t109 - t102 * t97) - t100 * t115, -g(1) * (-t97 * t108 + t109) - g(2) * (-t97 * t111 - t107) + t97 * t115, -t75, -g(1) * (t105 * t102 + t112) - g(2) * (t105 * t99 + t106) - g(3) * (t98 * pkin(2) - t101 * pkin(8) + pkin(6)) 0, 0, 0, 0, 0, 0, -g(1) * (t86 * t108 + t99 * t85) - g(2) * (-t102 * t85 + t86 * t111) - t86 * t115, -g(1) * (-t85 * t108 + t99 * t86) - g(2) * (-t102 * t86 - t85 * t111) + t85 * t115, -t75, -g(1) * (pkin(3) * t113 + t112) - g(2) * (-t99 * t110 + t84 * t111 + t91) - g(3) * (t101 * t103 + t98 * t84 + pkin(6)) + (-g(1) * (t101 * t84 - t110) - g(2) * (-pkin(7) - t89)) * t102, 0, 0, 0, 0, 0, 0, t73, t72, -t75, -g(1) * (t99 * t78 + t112) - g(2) * (t77 * t111 - t95 * t114 + t91) - g(3) * (t101 * t95 + t98 * t77 + pkin(6)) + (-g(1) * (t101 * t77 - t95 * t98) - g(2) * (-pkin(7) - t78)) * t102, 0, 0, 0, 0, 0, 0, t73, t72, -t75, -g(1) * (t99 * t76 + t112) - g(2) * (t74 * t111 - t87 * t114 + t91) - g(3) * (t101 * t87 + t98 * t74 + pkin(6)) + (-g(1) * (t101 * t74 - t87 * t98) - g(2) * (-pkin(7) - t76)) * t102;];
U_reg  = t1;
