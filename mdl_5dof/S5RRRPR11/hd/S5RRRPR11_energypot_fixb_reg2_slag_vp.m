% Calculate inertial parameters regressor of potential energy for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:13
% EndTime: 2019-12-31 21:35:13
% DurationCPUTime: 0.12s
% Computational Cost: add. (111->59), mult. (231->80), div. (0->0), fcn. (250->8), ass. (0->36)
t78 = sin(qJ(3));
t79 = sin(qJ(2));
t100 = t78 * t79;
t82 = cos(qJ(3));
t98 = t79 * t82;
t104 = pkin(3) * t98 + qJ(4) * t100;
t103 = g(3) * pkin(5);
t102 = g(3) * t79;
t101 = t79 * pkin(2) + pkin(5);
t80 = sin(qJ(1));
t99 = t79 * t80;
t84 = cos(qJ(1));
t97 = t79 * t84;
t83 = cos(qJ(2));
t96 = t80 * t83;
t95 = t84 * t78;
t94 = t84 * t82;
t93 = t84 * pkin(1) + t80 * pkin(6);
t92 = t80 * pkin(1) - t84 * pkin(6);
t91 = t84 * t83 * pkin(2) + pkin(7) * t97 + t93;
t90 = -t83 * pkin(7) + t101;
t89 = g(1) * t84 + g(2) * t80;
t88 = pkin(2) * t96 + pkin(7) * t99 + t92;
t62 = -t80 * t82 + t83 * t95;
t63 = t80 * t78 + t83 * t94;
t87 = t63 * pkin(3) + t62 * qJ(4) + t91;
t60 = t78 * t96 + t94;
t86 = g(1) * t62 + g(2) * t60 + g(3) * t100;
t61 = t82 * t96 - t95;
t85 = t61 * pkin(3) + t60 * qJ(4) + t88;
t81 = cos(qJ(5));
t77 = sin(qJ(5));
t64 = g(1) * t80 - g(2) * t84;
t57 = -g(3) * t83 + t89 * t79;
t56 = -g(1) * t63 - g(2) * t61 - g(3) * t98;
t1 = [0, 0, 0, 0, 0, 0, -t89, t64, -g(3), -t103, 0, 0, 0, 0, 0, 0, -t89 * t83 - t102, t57, -t64, -g(1) * t93 - g(2) * t92 - t103, 0, 0, 0, 0, 0, 0, t56, t86, -t57, -g(1) * t91 - g(2) * t88 - g(3) * t90, 0, 0, 0, 0, 0, 0, t56, -t57, -t86, -g(1) * t87 - g(2) * t85 - g(3) * (t90 + t104), 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t77 + t63 * t81) - g(2) * (t60 * t77 + t61 * t81) - (t77 * t78 + t81 * t82) * t102, -g(1) * (t62 * t81 - t63 * t77) - g(2) * (t60 * t81 - t61 * t77) - (-t77 * t82 + t78 * t81) * t102, t57, -g(1) * (t63 * pkin(4) - pkin(8) * t97 + t87) - g(2) * (t61 * pkin(4) - pkin(8) * t99 + t85) - g(3) * (pkin(4) * t98 + (-pkin(7) + pkin(8)) * t83 + t101 + t104);];
U_reg = t1;
