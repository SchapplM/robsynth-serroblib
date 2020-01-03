% Calculate inertial parameters regressor of potential energy for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:58
% EndTime: 2019-12-31 19:47:58
% DurationCPUTime: 0.12s
% Computational Cost: add. (96->65), mult. (161->78), div. (0->0), fcn. (156->8), ass. (0->39)
t73 = sin(qJ(1));
t72 = sin(qJ(2));
t83 = qJ(3) * t72;
t74 = cos(qJ(2));
t89 = t73 * t74;
t98 = pkin(2) * t89 + t73 * t83;
t97 = g(3) * pkin(5);
t69 = sin(pkin(8));
t96 = pkin(4) * t69;
t95 = g(3) * t74;
t94 = t72 * pkin(2) + pkin(5);
t68 = pkin(8) + qJ(5);
t61 = sin(t68);
t93 = t73 * t61;
t62 = cos(t68);
t92 = t73 * t62;
t91 = t73 * t69;
t70 = cos(pkin(8));
t90 = t73 * t70;
t75 = cos(qJ(1));
t88 = t75 * t61;
t87 = t75 * t62;
t86 = t75 * t69;
t85 = t75 * t70;
t84 = t75 * pkin(1) + t73 * pkin(6);
t82 = qJ(4) * t74;
t81 = t72 * t91;
t65 = t73 * pkin(1);
t80 = t65 + t98;
t79 = -t75 * pkin(6) + t65;
t78 = t84 + (pkin(2) * t74 + t83) * t75;
t77 = -t74 * qJ(3) + t94;
t76 = g(1) * t75 + g(2) * t73;
t71 = -pkin(7) - qJ(4);
t60 = t70 * pkin(4) + pkin(3);
t55 = g(1) * t73 - g(2) * t75;
t54 = g(3) * t72 + t76 * t74;
t53 = t76 * t72 - t95;
t1 = [0, 0, 0, 0, 0, 0, -t76, t55, -g(3), -t97, 0, 0, 0, 0, 0, 0, -t54, t53, -t55, -g(1) * t84 - g(2) * t79 - t97, 0, 0, 0, 0, 0, 0, -t55, t54, -t53, -g(1) * t78 - g(2) * (t79 + t98) - g(3) * t77, 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t86 + t90) - g(2) * (t81 - t85) + t69 * t95, -g(1) * (t72 * t85 - t91) - g(2) * (t72 * t90 + t86) + t70 * t95, -t54, -g(1) * (t73 * pkin(3) + t75 * t82 + t78) - g(2) * (t73 * t82 + (-pkin(3) - pkin(6)) * t75 + t80) - g(3) * (t72 * qJ(4) + t77), 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t88 + t92) - g(2) * (t72 * t93 - t87) + t61 * t95, -g(1) * (t72 * t87 - t93) - g(2) * (t72 * t92 + t88) + t62 * t95, -t54, -g(1) * (t73 * t60 + t78) - g(2) * (pkin(4) * t81 - t71 * t89 + t80) - g(3) * (-t72 * t71 + (-qJ(3) - t96) * t74 + t94) + (-g(1) * (-t71 * t74 + t72 * t96) - g(2) * (-pkin(6) - t60)) * t75;];
U_reg = t1;
