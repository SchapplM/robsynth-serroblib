% Calculate inertial parameters regressor of potential energy for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:29:50
% EndTime: 2019-03-09 04:29:50
% DurationCPUTime: 0.13s
% Computational Cost: add. (227->69), mult. (201->88), div. (0->0), fcn. (201->10), ass. (0->41)
t85 = qJ(2) + pkin(6);
t109 = g(3) * t85;
t87 = sin(qJ(3));
t108 = g(3) * t87;
t83 = qJ(1) + pkin(9);
t76 = sin(t83);
t86 = sin(qJ(4));
t107 = t76 * t86;
t90 = cos(qJ(3));
t106 = t76 * t90;
t78 = cos(t83);
t105 = t78 * t90;
t84 = -qJ(5) - pkin(8);
t104 = t84 * t87;
t103 = t86 * t90;
t89 = cos(qJ(4));
t102 = t89 * t90;
t88 = sin(qJ(1));
t101 = t88 * pkin(1) + t76 * pkin(2);
t91 = cos(qJ(1));
t100 = t91 * pkin(1) + t78 * pkin(2) + t76 * pkin(7);
t74 = t89 * pkin(4) + pkin(3);
t99 = t87 * t74 + t90 * t84 + t85;
t98 = -t78 * pkin(7) + t101;
t97 = pkin(3) * t90 + pkin(8) * t87;
t96 = g(1) * t78 + g(2) * t76;
t95 = -g(1) * t91 - g(2) * t88;
t82 = qJ(4) + pkin(10);
t75 = sin(t82);
t77 = cos(t82);
t59 = t75 * t106 + t78 * t77;
t61 = t75 * t105 - t76 * t77;
t94 = g(1) * t61 + g(2) * t59 + t75 * t108;
t93 = pkin(4) * t107 - t78 * t104 + t74 * t105 + t100;
t92 = -t76 * t104 + t74 * t106 + (-pkin(4) * t86 - pkin(7)) * t78 + t101;
t66 = g(1) * t76 - g(2) * t78;
t63 = -g(3) * t90 + t96 * t87;
t62 = t77 * t105 + t76 * t75;
t60 = t77 * t106 - t78 * t75;
t58 = -g(1) * t62 - g(2) * t60 - t77 * t108;
t1 = [0, 0, 0, 0, 0, 0, t95, g(1) * t88 - g(2) * t91, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t96, t66, -g(3), t95 * pkin(1) - t109, 0, 0, 0, 0, 0, 0, -t96 * t90 - t108, t63, -t66, -g(1) * t100 - g(2) * t98 - t109, 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t102 + t107) - g(2) * (t76 * t102 - t78 * t86) - t89 * t108, -g(1) * (-t78 * t103 + t76 * t89) - g(2) * (-t76 * t103 - t78 * t89) + t86 * t108, -t63, -g(1) * (t97 * t78 + t100) - g(2) * (t97 * t76 + t98) - g(3) * (t87 * pkin(3) - t90 * pkin(8) + t85) 0, 0, 0, 0, 0, 0, t58, t94, -t63, -g(1) * t93 - g(2) * t92 - g(3) * t99, 0, 0, 0, 0, 0, 0, t58, -t63, -t94, -g(1) * (t62 * pkin(5) + t61 * qJ(6) + t93) - g(2) * (t60 * pkin(5) + t59 * qJ(6) + t92) - g(3) * ((pkin(5) * t77 + qJ(6) * t75) * t87 + t99);];
U_reg  = t1;
