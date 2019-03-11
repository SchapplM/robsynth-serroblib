% Calculate inertial parameters regressor of potential energy for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:20
% EndTime: 2019-03-09 04:33:20
% DurationCPUTime: 0.13s
% Computational Cost: add. (209->57), mult. (222->69), div. (0->0), fcn. (228->8), ass. (0->37)
t83 = sin(qJ(3));
t86 = cos(qJ(3));
t108 = pkin(3) * t86 + pkin(8) * t83;
t85 = cos(qJ(4));
t101 = t83 * t85;
t82 = sin(qJ(4));
t103 = t82 * t83;
t107 = pkin(4) * t101 + qJ(5) * t103;
t81 = qJ(2) + pkin(6);
t104 = g(3) * t81;
t102 = t82 * t86;
t100 = t85 * t86;
t99 = qJ(6) * t83;
t98 = t83 * pkin(3) + t81;
t80 = qJ(1) + pkin(9);
t74 = sin(t80);
t75 = cos(t80);
t87 = cos(qJ(1));
t97 = t87 * pkin(1) + t75 * pkin(2) + t74 * pkin(7);
t84 = sin(qJ(1));
t96 = t84 * pkin(1) + t74 * pkin(2) - t75 * pkin(7);
t95 = t108 * t75 + t97;
t94 = g(1) * t75 + g(2) * t74;
t93 = -g(1) * t87 - g(2) * t84;
t92 = -t86 * pkin(8) + t98;
t91 = t108 * t74 + t96;
t59 = t74 * t102 + t75 * t85;
t61 = t75 * t102 - t74 * t85;
t90 = g(1) * t61 + g(2) * t59 + g(3) * t103;
t62 = t75 * t100 + t74 * t82;
t89 = t62 * pkin(4) + t61 * qJ(5) + t95;
t60 = t74 * t100 - t75 * t82;
t88 = t60 * pkin(4) + t59 * qJ(5) + t91;
t63 = g(1) * t74 - g(2) * t75;
t56 = -g(3) * t86 + t94 * t83;
t55 = -g(1) * t62 - g(2) * t60 - g(3) * t101;
t1 = [0, 0, 0, 0, 0, 0, t93, g(1) * t84 - g(2) * t87, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t94, t63, -g(3), t93 * pkin(1) - t104, 0, 0, 0, 0, 0, 0, -g(3) * t83 - t94 * t86, t56, -t63, -g(1) * t97 - g(2) * t96 - t104, 0, 0, 0, 0, 0, 0, t55, t90, -t56, -g(1) * t95 - g(2) * t91 - g(3) * t92, 0, 0, 0, 0, 0, 0, t55, -t56, -t90, -g(1) * t89 - g(2) * t88 - g(3) * (t92 + t107) 0, 0, 0, 0, 0, 0, t55, -t90, t56, -g(1) * (t62 * pkin(5) - t75 * t99 + t89) - g(2) * (t60 * pkin(5) - t74 * t99 + t88) - g(3) * (pkin(5) * t101 + (-pkin(8) + qJ(6)) * t86 + t98 + t107);];
U_reg  = t1;
