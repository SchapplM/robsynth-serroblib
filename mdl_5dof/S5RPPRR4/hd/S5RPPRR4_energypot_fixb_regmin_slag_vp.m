% Calculate minimal parameter regressor of potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:58
% EndTime: 2022-01-23 09:16:58
% DurationCPUTime: 0.08s
% Computational Cost: add. (73->42), mult. (91->64), div. (0->0), fcn. (97->10), ass. (0->31)
t81 = sin(pkin(8));
t100 = g(3) * t81;
t79 = pkin(9) + qJ(4);
t77 = qJ(5) + t79;
t73 = sin(t77);
t84 = sin(qJ(1));
t99 = t84 * t73;
t74 = cos(t77);
t98 = t84 * t74;
t75 = sin(t79);
t97 = t84 * t75;
t76 = cos(t79);
t96 = t84 * t76;
t80 = sin(pkin(9));
t95 = t84 * t80;
t82 = cos(pkin(9));
t94 = t84 * t82;
t85 = cos(qJ(1));
t93 = t85 * t73;
t92 = t85 * t74;
t91 = t85 * t75;
t90 = t85 * t76;
t89 = t85 * t80;
t88 = t85 * t82;
t87 = t85 * qJ(2);
t86 = -g(1) * t85 - g(2) * t84;
t83 = cos(pkin(8));
t78 = t84 * qJ(2);
t72 = g(1) * t84 - g(2) * t85;
t71 = pkin(2) * t83 + t81 * qJ(3) + pkin(1);
t1 = [0, t86, t72, t86 * t83 - t100, -t72, -g(1) * (t85 * pkin(1) + t78) - g(2) * (t84 * pkin(1) - t87) - g(3) * pkin(5), -g(1) * (t83 * t88 + t95) - g(2) * (t83 * t94 - t89) - t82 * t100, -g(1) * (-t83 * t89 + t94) - g(2) * (-t83 * t95 - t88) + t80 * t100, -g(1) * (t71 * t85 + t78) - g(2) * (t71 * t84 - t87) - g(3) * (t81 * pkin(2) - t83 * qJ(3) + pkin(5)), 0, 0, 0, 0, 0, -g(1) * (t83 * t90 + t97) - g(2) * (t83 * t96 - t91) - t76 * t100, -g(1) * (-t83 * t91 + t96) - g(2) * (-t83 * t97 - t90) + t75 * t100, 0, 0, 0, 0, 0, -g(1) * (t83 * t92 + t99) - g(2) * (t83 * t98 - t93) - t74 * t100, -g(1) * (-t83 * t93 + t98) - g(2) * (-t83 * t99 - t92) + t73 * t100;];
U_reg = t1;
