% Calculate inertial parameters regressor of potential energy for
% S5RPRRR6
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:39
% EndTime: 2019-12-31 19:01:39
% DurationCPUTime: 0.08s
% Computational Cost: add. (130->49), mult. (106->61), div. (0->0), fcn. (97->10), ass. (0->30)
t61 = qJ(3) + qJ(4);
t55 = sin(t61);
t79 = g(3) * t55;
t62 = qJ(2) + pkin(5);
t78 = g(3) * t62;
t56 = cos(t61);
t63 = sin(qJ(5));
t77 = t56 * t63;
t66 = cos(qJ(5));
t76 = t56 * t66;
t67 = cos(qJ(3));
t52 = t67 * pkin(3) + pkin(2);
t60 = qJ(1) + pkin(9);
t53 = sin(t60);
t54 = cos(t60);
t65 = sin(qJ(1));
t58 = t65 * pkin(1);
t69 = -pkin(7) - pkin(6);
t75 = t53 * t52 + t54 * t69 + t58;
t64 = sin(qJ(3));
t74 = t64 * pkin(3) + t62;
t68 = cos(qJ(1));
t59 = t68 * pkin(1);
t73 = t54 * t52 - t53 * t69 + t59;
t72 = pkin(4) * t56 + pkin(8) * t55;
t71 = g(1) * t54 + g(2) * t53;
t70 = -g(1) * t68 - g(2) * t65;
t47 = g(1) * t53 - g(2) * t54;
t46 = -g(3) * t56 + t71 * t55;
t1 = [0, 0, 0, 0, 0, 0, t70, g(1) * t65 - g(2) * t68, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t71, t47, -g(3), t70 * pkin(1) - t78, 0, 0, 0, 0, 0, 0, -g(3) * t64 - t71 * t67, -g(3) * t67 + t71 * t64, -t47, -g(1) * (t54 * pkin(2) + t53 * pkin(6) + t59) - g(2) * (t53 * pkin(2) - t54 * pkin(6) + t58) - t78, 0, 0, 0, 0, 0, 0, -t71 * t56 - t79, t46, -t47, -g(1) * t73 - g(2) * t75 - g(3) * t74, 0, 0, 0, 0, 0, 0, -g(1) * (t53 * t63 + t54 * t76) - g(2) * (t53 * t76 - t54 * t63) - t66 * t79, -g(1) * (t53 * t66 - t54 * t77) - g(2) * (-t53 * t77 - t54 * t66) + t63 * t79, -t46, -g(1) * (t72 * t54 + t73) - g(2) * (t72 * t53 + t75) - g(3) * (t55 * pkin(4) - t56 * pkin(8) + t74);];
U_reg = t1;
