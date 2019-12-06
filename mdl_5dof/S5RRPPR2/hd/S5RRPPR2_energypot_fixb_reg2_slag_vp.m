% Calculate inertial parameters regressor of potential energy for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:21
% EndTime: 2019-12-05 18:20:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (134->45), mult. (94->57), div. (0->0), fcn. (85->10), ass. (0->29)
t75 = pkin(6) + pkin(5);
t55 = qJ(3) + t75;
t74 = g(1) * t55;
t57 = sin(pkin(9));
t73 = g(1) * t57;
t56 = qJ(1) + qJ(2);
t51 = pkin(8) + t56;
t48 = sin(t51);
t72 = g(2) * t48;
t58 = cos(pkin(9));
t59 = sin(qJ(5));
t71 = t58 * t59;
t61 = cos(qJ(5));
t70 = t58 * t61;
t53 = cos(t56);
t62 = cos(qJ(1));
t69 = t62 * pkin(1) + pkin(2) * t53;
t49 = cos(t51);
t68 = t49 * pkin(3) + t48 * qJ(4) + t69;
t52 = sin(t56);
t60 = sin(qJ(1));
t67 = -t60 * pkin(1) - pkin(2) * t52;
t66 = pkin(4) * t58 + pkin(7) * t57;
t65 = -g(3) * t49 + t72;
t64 = g(2) * t60 - g(3) * t62;
t63 = t49 * qJ(4) + t67;
t44 = g(2) * t49 + g(3) * t48;
t43 = g(1) * t58 + t65 * t57;
t1 = [0, 0, 0, 0, 0, 0, t64, g(2) * t62 + g(3) * t60, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, g(2) * t52 - g(3) * t53, g(2) * t53 + g(3) * t52, -g(1), t64 * pkin(1) - g(1) * t75, 0, 0, 0, 0, 0, 0, t65, t44, -g(1), -g(2) * t67 - g(3) * t69 - t74, 0, 0, 0, 0, 0, 0, t65 * t58 - t73, -t43, -t44, -t74 - g(2) * (-t48 * pkin(3) + t63) - g(3) * t68, 0, 0, 0, 0, 0, 0, -t61 * t73 - g(2) * (-t48 * t70 + t49 * t59) - g(3) * (t48 * t59 + t49 * t70), t59 * t73 - g(2) * (t48 * t71 + t49 * t61) - g(3) * (t48 * t61 - t49 * t71), t43, -g(1) * (t57 * pkin(4) - t58 * pkin(7) + t55) - g(2) * t63 - g(3) * (t66 * t49 + t68) - (-pkin(3) - t66) * t72;];
U_reg = t1;
