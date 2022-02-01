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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:05:46
% EndTime: 2022-01-20 10:05:46
% DurationCPUTime: 0.10s
% Computational Cost: add. (134->45), mult. (94->57), div. (0->0), fcn. (85->10), ass. (0->28)
t76 = pkin(6) + pkin(5);
t57 = qJ(3) + t76;
t75 = g(3) * t57;
t59 = sin(pkin(9));
t74 = g(3) * t59;
t60 = cos(pkin(9));
t61 = sin(qJ(5));
t73 = t60 * t61;
t63 = cos(qJ(5));
t72 = t60 * t63;
t58 = qJ(1) + qJ(2);
t53 = sin(t58);
t62 = sin(qJ(1));
t71 = t62 * pkin(1) + pkin(2) * t53;
t54 = cos(t58);
t64 = cos(qJ(1));
t70 = t64 * pkin(1) + pkin(2) * t54;
t52 = pkin(8) + t58;
t48 = sin(t52);
t49 = cos(t52);
t69 = t49 * pkin(3) + t48 * qJ(4) + t70;
t68 = pkin(4) * t60 + pkin(7) * t59;
t67 = g(1) * t49 + g(2) * t48;
t66 = -g(1) * t64 - g(2) * t62;
t65 = t48 * pkin(3) - t49 * qJ(4) + t71;
t44 = g(1) * t48 - g(2) * t49;
t43 = -g(3) * t60 + t67 * t59;
t1 = [0, 0, 0, 0, 0, 0, t66, g(1) * t62 - g(2) * t64, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t54 - g(2) * t53, g(1) * t53 - g(2) * t54, -g(3), t66 * pkin(1) - g(3) * t76, 0, 0, 0, 0, 0, 0, -t67, t44, -g(3), -g(1) * t70 - g(2) * t71 - t75, 0, 0, 0, 0, 0, 0, -t67 * t60 - t74, t43, -t44, -g(1) * t69 - g(2) * t65 - t75, 0, 0, 0, 0, 0, 0, -g(1) * (t48 * t61 + t49 * t72) - g(2) * (t48 * t72 - t49 * t61) - t63 * t74, -g(1) * (t48 * t63 - t49 * t73) - g(2) * (-t48 * t73 - t49 * t63) + t61 * t74, -t43, -g(1) * (t68 * t49 + t69) - g(2) * (t68 * t48 + t65) - g(3) * (t59 * pkin(4) - t60 * pkin(7) + t57);];
U_reg = t1;
