% Calculate inertial parameters regressor of potential energy for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:36
% EndTime: 2019-12-31 20:30:36
% DurationCPUTime: 0.12s
% Computational Cost: add. (100->53), mult. (204->73), div. (0->0), fcn. (216->8), ass. (0->34)
t82 = sin(qJ(2));
t85 = cos(qJ(4));
t102 = t82 * t85;
t81 = sin(qJ(4));
t86 = cos(qJ(2));
t105 = t86 * t81 - t102;
t104 = g(3) * pkin(5);
t103 = g(3) * t105;
t83 = sin(qJ(1));
t101 = t83 * t86;
t87 = cos(qJ(1));
t99 = t86 * t87;
t98 = t87 * pkin(1) + t83 * pkin(6);
t97 = qJ(3) * t82;
t96 = t83 * pkin(1) - t87 * pkin(6);
t95 = pkin(2) * t99 + t87 * t97 + t98;
t94 = t82 * pkin(2) - t86 * qJ(3) + pkin(5);
t93 = g(1) * t87 + g(2) * t83;
t62 = t82 * t81 + t86 * t85;
t92 = pkin(2) * t101 + t83 * t97 + t96;
t91 = t82 * pkin(3) + t94;
t58 = t105 * t83;
t60 = -t87 * t102 + t81 * t99;
t90 = g(1) * t60 + g(2) * t58 + g(3) * t62;
t89 = pkin(3) * t101 + t87 * pkin(7) + t92;
t88 = pkin(3) * t99 - t83 * pkin(7) + t95;
t84 = cos(qJ(5));
t80 = sin(qJ(5));
t64 = g(1) * t83 - g(2) * t87;
t61 = t62 * t87;
t59 = t62 * t83;
t57 = -g(3) * t82 - t93 * t86;
t56 = -g(3) * t86 + t93 * t82;
t1 = [0, 0, 0, 0, 0, 0, -t93, t64, -g(3), -t104, 0, 0, 0, 0, 0, 0, t57, t56, -t64, -g(1) * t98 - g(2) * t96 - t104, 0, 0, 0, 0, 0, 0, t57, -t64, -t56, -g(1) * t95 - g(2) * t92 - g(3) * t94, 0, 0, 0, 0, 0, 0, -g(1) * t61 - g(2) * t59 + t103, t90, t64, -g(1) * t88 - g(2) * t89 - g(3) * t91, 0, 0, 0, 0, 0, 0, -g(1) * (t61 * t84 - t83 * t80) - g(2) * (t59 * t84 + t87 * t80) + t84 * t103, -g(1) * (-t61 * t80 - t83 * t84) - g(2) * (-t59 * t80 + t87 * t84) - t80 * t103, -t90, -g(1) * (t61 * pkin(4) + t60 * pkin(8) + t88) - g(2) * (t59 * pkin(4) + t58 * pkin(8) + t89) - g(3) * (-pkin(4) * t105 + t62 * pkin(8) + t91);];
U_reg = t1;
