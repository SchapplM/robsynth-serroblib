% Calculate inertial parameters regressor of potential energy for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:36
% EndTime: 2019-12-31 18:24:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (118->49), mult. (117->59), div. (0->0), fcn. (108->8), ass. (0->29)
t59 = qJ(1) + pkin(8);
t54 = sin(t59);
t62 = sin(qJ(3));
t73 = qJ(4) * t62;
t65 = cos(qJ(3));
t78 = t54 * t65;
t81 = pkin(3) * t78 + t54 * t73;
t60 = qJ(2) + pkin(5);
t80 = g(3) * t60;
t79 = g(3) * t65;
t55 = cos(t59);
t77 = t55 * t65;
t61 = sin(qJ(5));
t76 = t61 * t62;
t64 = cos(qJ(5));
t75 = t62 * t64;
t63 = sin(qJ(1));
t74 = t63 * pkin(1) + t54 * pkin(2);
t66 = cos(qJ(1));
t72 = t66 * pkin(1) + t55 * pkin(2) + t54 * pkin(6);
t71 = -t55 * pkin(6) + t74;
t70 = pkin(3) * t77 + t55 * t73 + t72;
t69 = g(1) * t55 + g(2) * t54;
t68 = -g(1) * t66 - g(2) * t63;
t67 = t62 * pkin(3) - t65 * qJ(4) + t60;
t45 = g(1) * t54 - g(2) * t55;
t44 = g(3) * t62 + t69 * t65;
t43 = t69 * t62 - t79;
t1 = [0, 0, 0, 0, 0, 0, t68, g(1) * t63 - g(2) * t66, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t69, t45, -g(3), t68 * pkin(1) - t80, 0, 0, 0, 0, 0, 0, -t44, t43, -t45, -g(1) * t72 - g(2) * t71 - t80, 0, 0, 0, 0, 0, 0, -t45, t44, -t43, -g(1) * t70 - g(2) * (t71 + t81) - g(3) * t67, 0, 0, 0, 0, 0, 0, -g(1) * (t54 * t64 + t55 * t76) - g(2) * (t54 * t76 - t55 * t64) + t61 * t79, -g(1) * (-t54 * t61 + t55 * t75) - g(2) * (t54 * t75 + t55 * t61) + t64 * t79, -t44, -g(1) * (t54 * pkin(4) + pkin(7) * t77 + t70) - g(2) * (pkin(7) * t78 + (-pkin(4) - pkin(6)) * t55 + t74 + t81) - g(3) * (t62 * pkin(7) + t67);];
U_reg = t1;
