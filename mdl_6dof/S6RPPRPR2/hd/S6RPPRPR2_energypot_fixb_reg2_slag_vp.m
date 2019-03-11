% Calculate inertial parameters regressor of potential energy for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:36
% EndTime: 2019-03-09 01:42:36
% DurationCPUTime: 0.15s
% Computational Cost: add. (191->62), mult. (147->70), div. (0->0), fcn. (135->10), ass. (0->39)
t79 = pkin(10) + qJ(4);
t74 = cos(t79);
t80 = qJ(1) + pkin(9);
t75 = cos(t80);
t100 = t74 * t75;
t72 = sin(t79);
t96 = qJ(5) * t72;
t106 = pkin(4) * t100 + t75 * t96;
t105 = g(3) * t74;
t83 = qJ(2) + pkin(6);
t104 = g(3) * t83;
t73 = sin(t80);
t103 = t73 * t74;
t85 = sin(qJ(6));
t102 = t73 * t85;
t87 = cos(qJ(6));
t101 = t73 * t87;
t99 = t75 * t85;
t98 = t75 * t87;
t82 = cos(pkin(10));
t71 = t82 * pkin(3) + pkin(2);
t88 = cos(qJ(1));
t78 = t88 * pkin(1);
t97 = t75 * t71 + t78;
t86 = sin(qJ(1));
t77 = t86 * pkin(1);
t84 = -pkin(7) - qJ(3);
t95 = t73 * t71 + t75 * t84 + t77;
t81 = sin(pkin(10));
t94 = t81 * pkin(3) + t83;
t93 = -t73 * t84 + t97;
t92 = pkin(4) * t103 + t73 * t96 + t95;
t91 = g(1) * t75 + g(2) * t73;
t90 = -g(1) * t88 - g(2) * t86;
t89 = t72 * pkin(4) - t74 * qJ(5) + t94;
t61 = g(1) * t73 - g(2) * t75;
t60 = g(3) * t72 + t91 * t74;
t59 = t91 * t72 - t105;
t1 = [0, 0, 0, 0, 0, 0, t90, g(1) * t86 - g(2) * t88, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t91, t61, -g(3), t90 * pkin(1) - t104, 0, 0, 0, 0, 0, 0, -g(3) * t81 - t91 * t82, -g(3) * t82 + t91 * t81, -t61, -g(1) * (t75 * pkin(2) + t73 * qJ(3) + t78) - g(2) * (t73 * pkin(2) - t75 * qJ(3) + t77) - t104, 0, 0, 0, 0, 0, 0, -t60, t59, -t61, -g(1) * t93 - g(2) * t95 - g(3) * t94, 0, 0, 0, 0, 0, 0, -t61, t60, -t59, -g(1) * (t93 + t106) - g(2) * t92 - g(3) * t89, 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t99 + t101) - g(2) * (t72 * t102 - t98) + t85 * t105, -g(1) * (t72 * t98 - t102) - g(2) * (t72 * t101 + t99) + t87 * t105, -t60, -g(1) * (pkin(8) * t100 + (pkin(5) - t84) * t73 + t97 + t106) - g(2) * (-t75 * pkin(5) + pkin(8) * t103 + t92) - g(3) * (t72 * pkin(8) + t89);];
U_reg  = t1;
