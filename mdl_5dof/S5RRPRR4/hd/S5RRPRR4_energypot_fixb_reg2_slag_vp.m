% Calculate inertial parameters regressor of potential energy for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:27
% EndTime: 2022-01-20 10:48:27
% DurationCPUTime: 0.07s
% Computational Cost: add. (118->42), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t70 = pkin(6) + pkin(5);
t57 = qJ(3) + t70;
t69 = g(3) * t57;
t59 = qJ(1) + qJ(2);
t52 = sin(t59);
t61 = sin(qJ(1));
t68 = t61 * pkin(1) + pkin(2) * t52;
t54 = cos(t59);
t63 = cos(qJ(1));
t67 = t63 * pkin(1) + pkin(2) * t54;
t50 = pkin(9) + t59;
t45 = sin(t50);
t46 = cos(t50);
t66 = g(1) * t46 + g(2) * t45;
t65 = -g(1) * t63 - g(2) * t61;
t64 = -pkin(8) - pkin(7);
t62 = cos(qJ(4));
t60 = sin(qJ(4));
t58 = qJ(4) + qJ(5);
t53 = cos(t58);
t51 = sin(t58);
t49 = t62 * pkin(4) + pkin(3);
t43 = g(1) * t45 - g(2) * t46;
t1 = [0, 0, 0, 0, 0, 0, t65, g(1) * t61 - g(2) * t63, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t54 - g(2) * t52, g(1) * t52 - g(2) * t54, -g(3), t65 * pkin(1) - g(3) * t70, 0, 0, 0, 0, 0, 0, -t66, t43, -g(3), -g(1) * t67 - g(2) * t68 - t69, 0, 0, 0, 0, 0, 0, -g(3) * t60 - t66 * t62, -g(3) * t62 + t66 * t60, -t43, -g(1) * (t46 * pkin(3) + t45 * pkin(7) + t67) - g(2) * (t45 * pkin(3) - t46 * pkin(7) + t68) - t69, 0, 0, 0, 0, 0, 0, -g(3) * t51 - t66 * t53, -g(3) * t53 + t66 * t51, -t43, -g(1) * (-t45 * t64 + t46 * t49 + t67) - g(2) * (t45 * t49 + t46 * t64 + t68) - g(3) * (t60 * pkin(4) + t57);];
U_reg = t1;
