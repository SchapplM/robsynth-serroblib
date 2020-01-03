% Calculate inertial parameters regressor of potential energy for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:09
% EndTime: 2019-12-31 17:19:09
% DurationCPUTime: 0.07s
% Computational Cost: add. (61->42), mult. (115->54), div. (0->0), fcn. (113->6), ass. (0->25)
t66 = g(3) * pkin(4);
t51 = sin(qJ(2));
t65 = g(3) * t51;
t49 = -qJ(4) - pkin(6);
t64 = t49 * t51;
t50 = sin(qJ(3));
t52 = sin(qJ(1));
t63 = t52 * t50;
t54 = cos(qJ(2));
t62 = t52 * t54;
t55 = cos(qJ(1));
t61 = t55 * t50;
t53 = cos(qJ(3));
t60 = t55 * t53;
t59 = t55 * pkin(1) + t52 * pkin(5);
t46 = t52 * pkin(1);
t58 = -t55 * pkin(5) + t46;
t57 = pkin(2) * t54 + pkin(6) * t51;
t56 = g(1) * t55 + g(2) * t52;
t44 = t53 * pkin(3) + pkin(2);
t43 = g(1) * t52 - g(2) * t55;
t42 = -g(3) * t54 + t56 * t51;
t41 = -g(1) * (t54 * t60 + t63) - g(2) * (t53 * t62 - t61) - t53 * t65;
t40 = -g(1) * (t52 * t53 - t54 * t61) - g(2) * (-t50 * t62 - t60) + t50 * t65;
t1 = [0, 0, 0, 0, 0, 0, -t56, t43, -g(3), -t66, 0, 0, 0, 0, 0, 0, -t56 * t54 - t65, t42, -t43, -g(1) * t59 - g(2) * t58 - t66, 0, 0, 0, 0, 0, 0, t41, t40, -t42, -g(1) * (t57 * t55 + t59) - g(2) * (t57 * t52 + t58) - g(3) * (t51 * pkin(2) - t54 * pkin(6) + pkin(4)), 0, 0, 0, 0, 0, 0, t41, t40, -t42, -g(1) * (pkin(3) * t63 + t59) - g(2) * (t44 * t62 - t52 * t64 + t46) - g(3) * (t51 * t44 + t54 * t49 + pkin(4)) + (-g(1) * (t44 * t54 - t64) - g(2) * (-pkin(3) * t50 - pkin(5))) * t55;];
U_reg = t1;
