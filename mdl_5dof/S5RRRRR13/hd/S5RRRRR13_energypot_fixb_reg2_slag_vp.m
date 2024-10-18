% Calculate inertial parameters regressor of potential energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:31:09
% EndTime: 2024-09-27 17:31:09
% DurationCPUTime: 0.03s
% Computational Cost: add. (177->61), mult. (117->83), div. (0->0), fcn. (106->16), ass. (0->38)
t92 = pkin(7) + pkin(6);
t76 = sin(pkin(5));
t91 = pkin(9) * t76;
t78 = sin(qJ(4));
t90 = t76 * t78;
t77 = cos(pkin(5));
t89 = t77 * t78;
t80 = cos(qJ(4));
t88 = t77 * t80;
t75 = qJ(1) + qJ(2);
t68 = sin(t75);
t79 = sin(qJ(1));
t87 = t79 * pkin(1) + pkin(2) * t68;
t70 = cos(t75);
t81 = cos(qJ(1));
t86 = t81 * pkin(1) + pkin(2) * t70;
t74 = qJ(4) + qJ(5);
t85 = pkin(8) + t92;
t71 = qJ(3) + t75;
t62 = sin(t71);
t63 = cos(t71);
t84 = g(1) * t62 - g(2) * t63;
t83 = -g(1) * t81 - g(2) * t79;
t82 = pkin(9) + pkin(10);
t69 = cos(t74);
t67 = sin(t74);
t66 = pkin(5) - t74;
t65 = pkin(5) + t74;
t64 = t80 * pkin(4) + pkin(3);
t59 = cos(t65);
t58 = sin(t66);
t57 = cos(t66) / 0.2e1;
t56 = sin(t65) / 0.2e1;
t55 = pkin(4) * t89 - t76 * t82;
t54 = t57 + t59 / 0.2e1;
t53 = t56 - t58 / 0.2e1;
t52 = -g(3) * t77 - t84 * t76;
t1 = [0, 0, 0, 0, 0, 0, t83, g(1) * t79 - g(2) * t81, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68, g(1) * t68 - g(2) * t70, -g(3), t83 * pkin(1) - g(3) * t92, 0, 0, 0, 0, 0, 0, -g(1) * t63 - g(2) * t62, t84, -g(3), -g(1) * t86 - g(2) * t87 - g(3) * t85, 0, 0, 0, 0, 0, 0, -g(1) * (-t62 * t89 + t63 * t80) - g(2) * (t62 * t80 + t63 * t89) - g(3) * t90, -g(1) * (-t62 * t88 - t63 * t78) - g(2) * (-t62 * t78 + t63 * t88) - g(3) * t76 * t80, t52, -g(1) * (t63 * pkin(3) + t62 * t91 + t86) - g(2) * (t62 * pkin(3) - t63 * t91 + t87) - g(3) * (t77 * pkin(9) + t85), 0, 0, 0, 0, 0, 0, -g(1) * (-t62 * t53 + t63 * t69) - g(2) * (t63 * t53 + t62 * t69) - g(3) * (t57 - t59 / 0.2e1), -g(1) * (-t62 * t54 - t63 * t67) - g(2) * (t63 * t54 - t62 * t67) - g(3) * (t56 + t58 / 0.2e1), t52, -g(1) * (-t62 * t55 + t63 * t64 + t86) - g(2) * (t63 * t55 + t62 * t64 + t87) - g(3) * (pkin(4) * t90 + t77 * t82 + t85);];
U_reg = t1;
