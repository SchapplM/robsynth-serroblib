% Calculate inertial parameters regressor of potential energy for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:30
% EndTime: 2019-12-31 20:01:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (129->49), mult. (156->61), div. (0->0), fcn. (155->8), ass. (0->34)
t75 = qJ(2) + pkin(8);
t71 = sin(t75);
t72 = cos(t75);
t99 = pkin(3) * t72 + pkin(7) * t71;
t98 = g(3) * pkin(5);
t95 = g(3) * t71;
t78 = sin(qJ(2));
t94 = t78 * pkin(2) + pkin(5);
t77 = sin(qJ(4));
t79 = sin(qJ(1));
t93 = t79 * t77;
t80 = cos(qJ(4));
t92 = t79 * t80;
t82 = cos(qJ(1));
t91 = t82 * t77;
t90 = t82 * t80;
t81 = cos(qJ(2));
t70 = t81 * pkin(2) + pkin(1);
t76 = -pkin(6) - qJ(3);
t89 = t79 * t70 + t82 * t76;
t88 = t82 * t70 - t79 * t76;
t87 = t99 * t79 + t89;
t86 = g(1) * t82 + g(2) * t79;
t85 = t71 * pkin(3) - t72 * pkin(7) + t94;
t84 = t99 * t82 + t88;
t55 = t72 * t93 + t90;
t57 = t72 * t91 - t92;
t83 = g(1) * t57 + g(2) * t55 + t77 * t95;
t59 = g(1) * t79 - g(2) * t82;
t58 = t72 * t90 + t93;
t56 = t72 * t92 - t91;
t54 = -g(3) * t72 + t86 * t71;
t53 = -g(1) * t58 - g(2) * t56 - t80 * t95;
t1 = [0, 0, 0, 0, 0, 0, -t86, t59, -g(3), -t98, 0, 0, 0, 0, 0, 0, -g(3) * t78 - t86 * t81, -g(3) * t81 + t86 * t78, -t59, -g(1) * (t82 * pkin(1) + t79 * pkin(6)) - g(2) * (t79 * pkin(1) - t82 * pkin(6)) - t98, 0, 0, 0, 0, 0, 0, -t86 * t72 - t95, t54, -t59, -g(1) * t88 - g(2) * t89 - g(3) * t94, 0, 0, 0, 0, 0, 0, t53, t83, -t54, -g(1) * t84 - g(2) * t87 - g(3) * t85, 0, 0, 0, 0, 0, 0, t53, -t54, -t83, -g(1) * (t58 * pkin(4) + t57 * qJ(5) + t84) - g(2) * (t56 * pkin(4) + t55 * qJ(5) + t87) - g(3) * ((pkin(4) * t80 + qJ(5) * t77) * t71 + t85);];
U_reg = t1;
