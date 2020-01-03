% Calculate inertial parameters regressor of potential energy for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:26
% EndTime: 2019-12-31 20:04:26
% DurationCPUTime: 0.13s
% Computational Cost: add. (88->49), mult. (168->58), div. (0->0), fcn. (167->6), ass. (0->31)
t72 = sin(qJ(1));
t71 = sin(qJ(2));
t84 = qJ(3) * t71;
t74 = cos(qJ(2));
t87 = t72 * t74;
t93 = pkin(2) * t87 + t72 * t84;
t75 = cos(qJ(1));
t78 = g(1) * t75 + g(2) * t72;
t92 = g(3) * pkin(5);
t89 = t71 * pkin(2) + pkin(5);
t70 = sin(qJ(4));
t88 = t71 * t70;
t86 = t74 * t75;
t85 = t75 * pkin(1) + t72 * pkin(6);
t83 = pkin(4) * t88;
t66 = t72 * pkin(1);
t82 = t66 + t93;
t81 = -t75 * pkin(6) + t66;
t80 = pkin(2) * t86 + t75 * t84 + t85;
t79 = -t74 * qJ(3) + t89;
t73 = cos(qJ(4));
t77 = t74 * t70 - t71 * t73;
t76 = t74 * t73 + t88;
t69 = -qJ(5) - pkin(7);
t63 = t73 * pkin(4) + pkin(3);
t58 = g(1) * t72 - g(2) * t75;
t57 = -g(3) * t71 - t74 * t78;
t56 = -g(3) * t74 + t71 * t78;
t55 = g(3) * t77 - t78 * t76;
t54 = g(3) * t76 + t78 * t77;
t1 = [0, 0, 0, 0, 0, 0, -t78, t58, -g(3), -t92, 0, 0, 0, 0, 0, 0, t57, t56, -t58, -g(1) * t85 - g(2) * t81 - t92, 0, 0, 0, 0, 0, 0, t57, -t58, -t56, -g(1) * t80 - g(2) * (t81 + t93) - g(3) * t79, 0, 0, 0, 0, 0, 0, t55, t54, t58, -g(1) * (pkin(3) * t86 - t72 * pkin(7) + t80) - g(2) * (pkin(3) * t87 + (-pkin(6) + pkin(7)) * t75 + t82) - g(3) * (t71 * pkin(3) + t79), 0, 0, 0, 0, 0, 0, t55, t54, t58, -g(1) * (t72 * t69 + t80) - g(2) * (t63 * t87 + t72 * t83 + t82) - g(3) * (t71 * t63 + (-pkin(4) * t70 - qJ(3)) * t74 + t89) + (-g(1) * (t63 * t74 + t83) - g(2) * (-pkin(6) - t69)) * t75;];
U_reg = t1;
