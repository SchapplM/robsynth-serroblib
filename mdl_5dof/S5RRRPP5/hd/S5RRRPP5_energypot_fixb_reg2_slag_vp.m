% Calculate inertial parameters regressor of potential energy for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:27
% EndTime: 2019-12-31 20:58:27
% DurationCPUTime: 0.07s
% Computational Cost: add. (105->43), mult. (117->48), div. (0->0), fcn. (104->6), ass. (0->25)
t70 = cos(qJ(1));
t66 = qJ(2) + qJ(3);
t62 = sin(t66);
t76 = qJ(4) * t62;
t63 = cos(t66);
t78 = t63 * t70;
t82 = pkin(3) * t78 + t70 * t76;
t81 = g(3) * pkin(5);
t67 = sin(qJ(2));
t80 = t67 * pkin(2) + pkin(5);
t68 = sin(qJ(1));
t79 = t63 * t68;
t69 = cos(qJ(2));
t60 = t69 * pkin(2) + pkin(1);
t71 = -pkin(7) - pkin(6);
t77 = t68 * t60 + t70 * t71;
t55 = t70 * t60;
t75 = -t68 * t71 + t55;
t74 = pkin(3) * t79 + t68 * t76 + t77;
t73 = g(1) * t70 + g(2) * t68;
t72 = t62 * pkin(3) - t63 * qJ(4) + t80;
t51 = g(1) * t68 - g(2) * t70;
t50 = -g(3) * t62 - t73 * t63;
t49 = -g(3) * t63 + t73 * t62;
t1 = [0, 0, 0, 0, 0, 0, -t73, t51, -g(3), -t81, 0, 0, 0, 0, 0, 0, -g(3) * t67 - t73 * t69, -g(3) * t69 + t73 * t67, -t51, -g(1) * (t70 * pkin(1) + t68 * pkin(6)) - g(2) * (t68 * pkin(1) - t70 * pkin(6)) - t81, 0, 0, 0, 0, 0, 0, t50, t49, -t51, -g(1) * t75 - g(2) * t77 - g(3) * t80, 0, 0, 0, 0, 0, 0, t50, -t51, -t49, -g(1) * (t75 + t82) - g(2) * t74 - g(3) * t72, 0, 0, 0, 0, 0, 0, t50, -t49, t51, -g(1) * (pkin(4) * t78 + t55 + (-qJ(5) - t71) * t68 + t82) - g(2) * (pkin(4) * t79 + t70 * qJ(5) + t74) - g(3) * (t62 * pkin(4) + t72);];
U_reg = t1;
