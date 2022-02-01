% Calculate inertial parameters regressor of potential energy for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:03
% EndTime: 2022-01-23 09:19:03
% DurationCPUTime: 0.07s
% Computational Cost: add. (118->42), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t69 = qJ(2) + pkin(5);
t57 = pkin(6) + t69;
t70 = g(3) * t57;
t59 = qJ(1) + pkin(8);
t51 = sin(t59);
t63 = sin(qJ(1));
t68 = t63 * pkin(1) + pkin(2) * t51;
t53 = cos(t59);
t64 = cos(qJ(1));
t67 = t64 * pkin(1) + pkin(2) * t53;
t54 = qJ(3) + t59;
t47 = sin(t54);
t48 = cos(t54);
t66 = g(1) * t48 + g(2) * t47;
t65 = -g(1) * t64 - g(2) * t63;
t62 = -pkin(7) - qJ(4);
t61 = cos(pkin(9));
t60 = sin(pkin(9));
t58 = pkin(9) + qJ(5);
t52 = cos(t58);
t50 = sin(t58);
t49 = t61 * pkin(4) + pkin(3);
t43 = g(1) * t47 - g(2) * t48;
t1 = [0, 0, 0, 0, 0, 0, t65, g(1) * t63 - g(2) * t64, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t53 - g(2) * t51, g(1) * t51 - g(2) * t53, -g(3), t65 * pkin(1) - g(3) * t69, 0, 0, 0, 0, 0, 0, -t66, t43, -g(3), -g(1) * t67 - g(2) * t68 - t70, 0, 0, 0, 0, 0, 0, -g(3) * t60 - t66 * t61, -g(3) * t61 + t66 * t60, -t43, -g(1) * (t48 * pkin(3) + t47 * qJ(4) + t67) - g(2) * (t47 * pkin(3) - t48 * qJ(4) + t68) - t70, 0, 0, 0, 0, 0, 0, -g(3) * t50 - t66 * t52, -g(3) * t52 + t66 * t50, -t43, -g(1) * (-t47 * t62 + t48 * t49 + t67) - g(2) * (t47 * t49 + t48 * t62 + t68) - g(3) * (t60 * pkin(4) + t57);];
U_reg = t1;
