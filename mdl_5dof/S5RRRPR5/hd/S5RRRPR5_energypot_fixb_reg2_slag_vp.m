% Calculate inertial parameters regressor of potential energy for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:16
% EndTime: 2019-12-31 21:14:16
% DurationCPUTime: 0.10s
% Computational Cost: add. (130->50), mult. (119->64), div. (0->0), fcn. (110->10), ass. (0->31)
t92 = g(3) * pkin(5);
t80 = -pkin(7) - pkin(6);
t73 = qJ(2) + qJ(3);
t66 = pkin(9) + t73;
t62 = sin(t66);
t91 = g(3) * t62;
t75 = sin(qJ(2));
t90 = t75 * pkin(2) + pkin(5);
t78 = cos(qJ(2));
t65 = t78 * pkin(2) + pkin(1);
t74 = sin(qJ(5));
t76 = sin(qJ(1));
t89 = t76 * t74;
t77 = cos(qJ(5));
t88 = t76 * t77;
t79 = cos(qJ(1));
t87 = t79 * t74;
t86 = t79 * t77;
t68 = cos(t73);
t59 = pkin(3) * t68 + t65;
t72 = -qJ(4) + t80;
t85 = t76 * t59 + t79 * t72;
t67 = sin(t73);
t84 = pkin(3) * t67 + t90;
t83 = t79 * t59 - t76 * t72;
t63 = cos(t66);
t82 = pkin(4) * t63 + pkin(8) * t62;
t81 = g(1) * t79 + g(2) * t76;
t60 = g(1) * t76 - g(2) * t79;
t56 = -g(3) * t63 + t81 * t62;
t1 = [0, 0, 0, 0, 0, 0, -t81, t60, -g(3), -t92, 0, 0, 0, 0, 0, 0, -g(3) * t75 - t81 * t78, -g(3) * t78 + t81 * t75, -t60, -g(1) * (t79 * pkin(1) + t76 * pkin(6)) - g(2) * (t76 * pkin(1) - t79 * pkin(6)) - t92, 0, 0, 0, 0, 0, 0, -g(3) * t67 - t81 * t68, -g(3) * t68 + t81 * t67, -t60, -g(1) * (t79 * t65 - t76 * t80) - g(2) * (t76 * t65 + t79 * t80) - g(3) * t90, 0, 0, 0, 0, 0, 0, -t81 * t63 - t91, t56, -t60, -g(1) * t83 - g(2) * t85 - g(3) * t84, 0, 0, 0, 0, 0, 0, -g(1) * (t63 * t86 + t89) - g(2) * (t63 * t88 - t87) - t77 * t91, -g(1) * (-t63 * t87 + t88) - g(2) * (-t63 * t89 - t86) + t74 * t91, -t56, -g(1) * (t82 * t79 + t83) - g(2) * (t82 * t76 + t85) - g(3) * (t62 * pkin(4) - t63 * pkin(8) + t84);];
U_reg = t1;
