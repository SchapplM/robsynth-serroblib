% Calculate inertial parameters regressor of potential energy for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:36
% EndTime: 2019-12-31 18:45:36
% DurationCPUTime: 0.12s
% Computational Cost: add. (136->46), mult. (143->58), div. (0->0), fcn. (142->8), ass. (0->31)
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t87 = pkin(3) * t71 + pkin(7) * t68;
t66 = qJ(2) + pkin(5);
t84 = g(3) * t66;
t83 = g(3) * t68;
t67 = sin(qJ(4));
t82 = t67 * t71;
t70 = cos(qJ(4));
t81 = t70 * t71;
t65 = qJ(1) + pkin(8);
t59 = sin(t65);
t60 = cos(t65);
t72 = cos(qJ(1));
t80 = t72 * pkin(1) + t60 * pkin(2) + t59 * pkin(6);
t69 = sin(qJ(1));
t79 = t69 * pkin(1) + t59 * pkin(2) - t60 * pkin(6);
t78 = t87 * t60 + t80;
t77 = g(1) * t60 + g(2) * t59;
t76 = -g(1) * t72 - g(2) * t69;
t75 = t68 * pkin(3) - t71 * pkin(7) + t66;
t74 = t87 * t59 + t79;
t46 = t59 * t82 + t60 * t70;
t48 = -t59 * t70 + t60 * t82;
t73 = g(1) * t48 + g(2) * t46 + t67 * t83;
t50 = g(1) * t59 - g(2) * t60;
t49 = t59 * t67 + t60 * t81;
t47 = t59 * t81 - t60 * t67;
t45 = -g(3) * t71 + t77 * t68;
t44 = -g(1) * t49 - g(2) * t47 - t70 * t83;
t1 = [0, 0, 0, 0, 0, 0, t76, g(1) * t69 - g(2) * t72, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t77, t50, -g(3), pkin(1) * t76 - t84, 0, 0, 0, 0, 0, 0, -t71 * t77 - t83, t45, -t50, -g(1) * t80 - g(2) * t79 - t84, 0, 0, 0, 0, 0, 0, t44, t73, -t45, -g(1) * t78 - g(2) * t74 - g(3) * t75, 0, 0, 0, 0, 0, 0, t44, -t45, -t73, -g(1) * (t49 * pkin(4) + t48 * qJ(5) + t78) - g(2) * (t47 * pkin(4) + t46 * qJ(5) + t74) - g(3) * ((pkin(4) * t70 + qJ(5) * t67) * t68 + t75);];
U_reg = t1;
