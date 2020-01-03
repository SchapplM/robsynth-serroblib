% Calculate inertial parameters regressor of potential energy for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:57
% EndTime: 2019-12-31 19:29:57
% DurationCPUTime: 0.10s
% Computational Cost: add. (117->46), mult. (133->56), div. (0->0), fcn. (126->8), ass. (0->29)
t80 = cos(qJ(1));
t73 = qJ(2) + pkin(8);
t69 = sin(t73);
t87 = qJ(4) * t69;
t70 = cos(t73);
t89 = t70 * t80;
t95 = pkin(3) * t89 + t80 * t87;
t77 = sin(qJ(1));
t84 = g(1) * t80 + g(2) * t77;
t94 = g(3) * pkin(5);
t76 = sin(qJ(2));
t91 = t76 * pkin(2) + pkin(5);
t90 = t70 * t77;
t79 = cos(qJ(2));
t68 = t79 * pkin(2) + pkin(1);
t74 = -pkin(6) - qJ(3);
t88 = t77 * t68 + t80 * t74;
t65 = t80 * t68;
t86 = -t77 * t74 + t65;
t85 = pkin(3) * t90 + t77 * t87 + t88;
t75 = sin(qJ(5));
t78 = cos(qJ(5));
t83 = t69 * t78 - t70 * t75;
t82 = t69 * t75 + t70 * t78;
t81 = t69 * pkin(3) - t70 * qJ(4) + t91;
t61 = g(1) * t77 - g(2) * t80;
t58 = -g(3) * t69 - t84 * t70;
t57 = -g(3) * t70 + t84 * t69;
t1 = [0, 0, 0, 0, 0, 0, -t84, t61, -g(3), -t94, 0, 0, 0, 0, 0, 0, -g(3) * t76 - t84 * t79, -g(3) * t79 + t84 * t76, -t61, -g(1) * (t80 * pkin(1) + t77 * pkin(6)) - g(2) * (t77 * pkin(1) - t80 * pkin(6)) - t94, 0, 0, 0, 0, 0, 0, t58, t57, -t61, -g(1) * t86 - g(2) * t88 - g(3) * t91, 0, 0, 0, 0, 0, 0, t58, -t61, -t57, -g(1) * (t86 + t95) - g(2) * t85 - g(3) * t81, 0, 0, 0, 0, 0, 0, -g(3) * t83 - t84 * t82, g(3) * t82 - t84 * t83, t61, -g(1) * (pkin(4) * t89 + t65 + (-pkin(7) - t74) * t77 + t95) - g(2) * (pkin(4) * t90 + t80 * pkin(7) + t85) - g(3) * (t69 * pkin(4) + t81);];
U_reg = t1;
