% Calculate inertial parameters regressor of potential energy for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:12
% EndTime: 2019-12-05 18:12:12
% DurationCPUTime: 0.10s
% Computational Cost: add. (118->45), mult. (99->54), div. (0->0), fcn. (86->10), ass. (0->25)
t81 = g(3) * pkin(5);
t73 = sin(pkin(9));
t80 = t73 * pkin(2) + pkin(5);
t74 = cos(pkin(9));
t62 = t74 * pkin(2) + pkin(1);
t75 = -pkin(6) - qJ(2);
t72 = pkin(9) + qJ(3);
t71 = -pkin(7) + t75;
t64 = sin(t72);
t79 = pkin(3) * t64 + t80;
t65 = cos(t72);
t54 = pkin(3) * t65 + t62;
t66 = qJ(4) + t72;
t76 = sin(qJ(1));
t77 = cos(qJ(1));
t78 = g(1) * t77 + g(2) * t76;
t67 = -pkin(8) + t71;
t63 = qJ(5) + t66;
t61 = cos(t66);
t60 = sin(t66);
t57 = cos(t63);
t56 = sin(t63);
t55 = g(1) * t76 - g(2) * t77;
t53 = pkin(4) * t61 + t54;
t1 = [0, 0, 0, 0, 0, 0, -t78, t55, -g(3), -t81, 0, 0, 0, 0, 0, 0, -g(3) * t73 - t78 * t74, -g(3) * t74 + t78 * t73, -t55, -g(1) * (t77 * pkin(1) + t76 * qJ(2)) - g(2) * (t76 * pkin(1) - t77 * qJ(2)) - t81, 0, 0, 0, 0, 0, 0, -g(3) * t64 - t78 * t65, -g(3) * t65 + t78 * t64, -t55, -g(1) * (t77 * t62 - t76 * t75) - g(2) * (t76 * t62 + t77 * t75) - g(3) * t80, 0, 0, 0, 0, 0, 0, -g(3) * t60 - t78 * t61, -g(3) * t61 + t78 * t60, -t55, -g(1) * (t77 * t54 - t76 * t71) - g(2) * (t76 * t54 + t77 * t71) - g(3) * t79, 0, 0, 0, 0, 0, 0, -g(3) * t56 - t78 * t57, -g(3) * t57 + t78 * t56, -t55, -g(1) * (t77 * t53 - t76 * t67) - g(2) * (t76 * t53 + t77 * t67) - g(3) * (pkin(4) * t60 + t79);];
U_reg = t1;
