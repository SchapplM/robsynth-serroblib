% Calculate inertial parameters regressor of potential energy for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR14_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:15
% EndTime: 2019-12-31 18:35:15
% DurationCPUTime: 0.09s
% Computational Cost: add. (89->50), mult. (112->58), div. (0->0), fcn. (103->8), ass. (0->30)
t53 = qJ(3) + pkin(8);
t47 = sin(t53);
t48 = cos(t53);
t77 = pkin(4) * t47 - pkin(7) * t48;
t76 = g(3) * pkin(5);
t75 = pkin(2) + pkin(5);
t56 = sin(qJ(3));
t74 = pkin(3) * t56;
t71 = g(3) * t48;
t55 = sin(qJ(5));
t57 = sin(qJ(1));
t70 = t57 * t55;
t58 = cos(qJ(5));
t69 = t57 * t58;
t60 = cos(qJ(1));
t68 = t60 * t55;
t67 = t60 * t58;
t66 = t60 * pkin(1) + t57 * qJ(2);
t59 = cos(qJ(3));
t65 = t59 * pkin(3) + t75;
t64 = t57 * t74 + t66;
t63 = -qJ(2) - t74;
t50 = t57 * pkin(1);
t54 = -qJ(4) - pkin(6);
t62 = -t57 * t54 + t50;
t61 = -t60 * qJ(2) + t50;
t44 = g(1) * t57 - g(2) * t60;
t45 = g(1) * t60 + g(2) * t57;
t43 = -g(3) * t47 + t44 * t48;
t1 = [0, 0, 0, 0, 0, 0, -t45, t44, -g(3), -t76, 0, 0, 0, 0, 0, 0, -g(3), t45, -t44, -g(1) * t66 - g(2) * t61 - t76, 0, 0, 0, 0, 0, 0, -g(3) * t59 - t44 * t56, g(3) * t56 - t44 * t59, -t45, -g(1) * (t60 * pkin(6) + t66) - g(2) * (t57 * pkin(6) + t61) - g(3) * t75, 0, 0, 0, 0, 0, 0, -t44 * t47 - t71, -t43, -t45, -g(1) * (-t60 * t54 + t64) - g(2) * (t63 * t60 + t62) - g(3) * t65, 0, 0, 0, 0, 0, 0, -g(1) * (t47 * t69 + t68) - g(2) * (-t47 * t67 + t70) - t58 * t71, -g(1) * (-t47 * t70 + t67) - g(2) * (t47 * t68 + t69) + t55 * t71, t43, -g(1) * (t77 * t57 + t64) - g(2) * t62 - g(3) * (t48 * pkin(4) + t47 * pkin(7) + t65) + (g(1) * t54 - g(2) * (t63 - t77)) * t60;];
U_reg = t1;
