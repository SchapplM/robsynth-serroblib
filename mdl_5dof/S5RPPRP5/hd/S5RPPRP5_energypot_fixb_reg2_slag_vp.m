% Calculate inertial parameters regressor of potential energy for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:42
% EndTime: 2019-12-31 17:53:42
% DurationCPUTime: 0.08s
% Computational Cost: add. (92->46), mult. (183->60), div. (0->0), fcn. (188->6), ass. (0->32)
t92 = g(3) * pkin(5);
t72 = sin(pkin(7));
t76 = cos(qJ(4));
t91 = t72 * t76;
t73 = cos(pkin(7));
t75 = sin(qJ(1));
t90 = t73 * t75;
t77 = cos(qJ(1));
t89 = t73 * t77;
t88 = t77 * pkin(1) + t75 * qJ(2);
t87 = qJ(3) * t72;
t86 = t75 * pkin(1) - t77 * qJ(2);
t85 = pkin(2) * t89 + t77 * t87 + t88;
t84 = t72 * pkin(2) - t73 * qJ(3) + pkin(5);
t83 = g(1) * t77 + g(2) * t75;
t74 = sin(qJ(4));
t54 = t72 * t74 + t73 * t76;
t82 = t72 * pkin(3) + t84;
t81 = pkin(2) * t90 + t75 * t87 + t86;
t50 = t74 * t90 - t75 * t91;
t52 = t74 * t89 - t77 * t91;
t80 = g(1) * t52 + g(2) * t50 + g(3) * t54;
t79 = pkin(3) * t89 - t75 * pkin(6) + t85;
t78 = pkin(3) * t90 + t77 * pkin(6) + t81;
t56 = g(1) * t75 - g(2) * t77;
t55 = -t73 * t74 + t91;
t53 = t54 * t77;
t51 = t54 * t75;
t49 = -g(3) * t72 - t83 * t73;
t48 = -g(3) * t73 + t83 * t72;
t47 = -g(1) * t53 - g(2) * t51 - g(3) * t55;
t1 = [0, 0, 0, 0, 0, 0, -t83, t56, -g(3), -t92, 0, 0, 0, 0, 0, 0, t49, t48, -t56, -g(1) * t88 - g(2) * t86 - t92, 0, 0, 0, 0, 0, 0, t49, -t56, -t48, -g(1) * t85 - g(2) * t81 - g(3) * t84, 0, 0, 0, 0, 0, 0, t47, t80, t56, -g(1) * t79 - g(2) * t78 - g(3) * t82, 0, 0, 0, 0, 0, 0, t47, t56, -t80, -g(1) * (t53 * pkin(4) + t52 * qJ(5) + t79) - g(2) * (t51 * pkin(4) + t50 * qJ(5) + t78) - g(3) * (t55 * pkin(4) + t54 * qJ(5) + t82);];
U_reg = t1;
