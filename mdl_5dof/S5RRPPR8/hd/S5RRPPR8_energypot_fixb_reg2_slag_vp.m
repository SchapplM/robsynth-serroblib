% Calculate inertial parameters regressor of potential energy for
% S5RRPPR8
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:14
% EndTime: 2019-12-31 19:39:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (100->54), mult. (168->66), div. (0->0), fcn. (167->8), ass. (0->34)
t76 = sin(qJ(1));
t75 = sin(qJ(2));
t89 = qJ(3) * t75;
t77 = cos(qJ(2));
t92 = t76 * t77;
t98 = pkin(2) * t92 + t76 * t89;
t78 = cos(qJ(1));
t83 = g(1) * t78 + g(2) * t76;
t97 = g(3) * pkin(5);
t94 = t75 * pkin(2) + pkin(5);
t72 = sin(pkin(8));
t93 = t75 * t72;
t91 = t77 * t78;
t90 = t78 * pkin(1) + t76 * pkin(6);
t88 = pkin(4) * t93;
t68 = t76 * pkin(1);
t87 = t68 + t98;
t86 = -t78 * pkin(6) + t68;
t85 = pkin(2) * t91 + t78 * t89 + t90;
t84 = -t77 * qJ(3) + t94;
t71 = pkin(8) + qJ(5);
t64 = sin(t71);
t65 = cos(t71);
t82 = t77 * t64 - t75 * t65;
t81 = t75 * t64 + t77 * t65;
t73 = cos(pkin(8));
t80 = t77 * t72 - t75 * t73;
t79 = t77 * t73 + t93;
t74 = -pkin(7) - qJ(4);
t63 = t73 * pkin(4) + pkin(3);
t58 = g(1) * t76 - g(2) * t78;
t57 = -g(3) * t75 - t83 * t77;
t56 = -g(3) * t77 + t83 * t75;
t1 = [0, 0, 0, 0, 0, 0, -t83, t58, -g(3), -t97, 0, 0, 0, 0, 0, 0, t57, t56, -t58, -g(1) * t90 - g(2) * t86 - t97, 0, 0, 0, 0, 0, 0, t57, -t58, -t56, -g(1) * t85 - g(2) * (t86 + t98) - g(3) * t84, 0, 0, 0, 0, 0, 0, g(3) * t80 - t83 * t79, g(3) * t79 + t83 * t80, t58, -g(1) * (pkin(3) * t91 - t76 * qJ(4) + t85) - g(2) * (pkin(3) * t92 + (-pkin(6) + qJ(4)) * t78 + t87) - g(3) * (t75 * pkin(3) + t84), 0, 0, 0, 0, 0, 0, g(3) * t82 - t83 * t81, g(3) * t81 + t83 * t82, t58, -g(1) * (t76 * t74 + t85) - g(2) * (t63 * t92 + t76 * t88 + t87) - g(3) * (t75 * t63 + (-pkin(4) * t72 - qJ(3)) * t77 + t94) + (-g(1) * (t63 * t77 + t88) - g(2) * (-pkin(6) - t74)) * t78;];
U_reg = t1;
