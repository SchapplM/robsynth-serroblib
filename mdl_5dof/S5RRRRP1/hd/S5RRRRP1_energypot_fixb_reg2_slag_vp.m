% Calculate inertial parameters regressor of potential energy for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:45:59
% EndTime: 2019-12-05 18:45:59
% DurationCPUTime: 0.07s
% Computational Cost: add. (112->42), mult. (99->50), div. (0->0), fcn. (86->8), ass. (0->24)
t81 = g(3) * pkin(5);
t77 = -pkin(7) - pkin(6);
t73 = sin(qJ(2));
t80 = t73 * pkin(2) + pkin(5);
t75 = cos(qJ(2));
t63 = t75 * pkin(2) + pkin(1);
t72 = qJ(2) + qJ(3);
t71 = -pkin(8) + t77;
t64 = sin(t72);
t79 = pkin(3) * t64 + t80;
t65 = cos(t72);
t57 = pkin(3) * t65 + t63;
t74 = sin(qJ(1));
t76 = cos(qJ(1));
t78 = g(1) * t76 + g(2) * t74;
t67 = qJ(4) + t72;
t66 = -qJ(5) + t71;
t62 = cos(t67);
t61 = sin(t67);
t58 = g(1) * t74 - g(2) * t76;
t56 = pkin(4) * t62 + t57;
t55 = -g(3) * t61 - t78 * t62;
t54 = -g(3) * t62 + t78 * t61;
t1 = [0, 0, 0, 0, 0, 0, -t78, t58, -g(3), -t81, 0, 0, 0, 0, 0, 0, -g(3) * t73 - t78 * t75, -g(3) * t75 + t78 * t73, -t58, -g(1) * (t76 * pkin(1) + t74 * pkin(6)) - g(2) * (t74 * pkin(1) - t76 * pkin(6)) - t81, 0, 0, 0, 0, 0, 0, -g(3) * t64 - t78 * t65, -g(3) * t65 + t78 * t64, -t58, -g(1) * (t76 * t63 - t74 * t77) - g(2) * (t74 * t63 + t76 * t77) - g(3) * t80, 0, 0, 0, 0, 0, 0, t55, t54, -t58, -g(1) * (t76 * t57 - t74 * t71) - g(2) * (t74 * t57 + t76 * t71) - g(3) * t79, 0, 0, 0, 0, 0, 0, t55, t54, -t58, -g(1) * (t76 * t56 - t74 * t66) - g(2) * (t74 * t56 + t76 * t66) - g(3) * (pkin(4) * t61 + t79);];
U_reg = t1;
