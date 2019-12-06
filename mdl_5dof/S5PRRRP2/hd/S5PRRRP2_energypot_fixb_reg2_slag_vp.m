% Calculate inertial parameters regressor of potential energy for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:56
% EndTime: 2019-12-05 16:41:56
% DurationCPUTime: 0.06s
% Computational Cost: add. (121->38), mult. (81->41), div. (0->0), fcn. (68->8), ass. (0->24)
t61 = pkin(5) + qJ(1);
t48 = pkin(6) + t61;
t62 = g(3) * t48;
t49 = pkin(8) + qJ(2);
t43 = sin(t49);
t50 = sin(pkin(8));
t60 = t50 * pkin(1) + pkin(2) * t43;
t44 = cos(t49);
t51 = cos(pkin(8));
t59 = t51 * pkin(1) + pkin(2) * t44;
t45 = qJ(3) + t49;
t41 = sin(t45);
t42 = cos(t45);
t58 = t42 * pkin(3) + t41 * pkin(7) + t59;
t57 = g(1) * t42 + g(2) * t41;
t56 = -g(1) * t51 - g(2) * t50;
t52 = sin(qJ(4));
t53 = cos(qJ(4));
t55 = pkin(4) * t53 + qJ(5) * t52;
t54 = t41 * pkin(3) - t42 * pkin(7) + t60;
t34 = g(1) * t41 - g(2) * t42;
t33 = -g(3) * t52 - t57 * t53;
t32 = -g(3) * t53 + t57 * t52;
t1 = [0, 0, 0, 0, 0, 0, t56, g(1) * t50 - g(2) * t51, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t44 - g(2) * t43, g(1) * t43 - g(2) * t44, -g(3), t56 * pkin(1) - g(3) * t61, 0, 0, 0, 0, 0, 0, -t57, t34, -g(3), -g(1) * t59 - g(2) * t60 - t62, 0, 0, 0, 0, 0, 0, t33, t32, -t34, -g(1) * t58 - g(2) * t54 - t62, 0, 0, 0, 0, 0, 0, t33, -t34, -t32, -g(1) * (t55 * t42 + t58) - g(2) * (t55 * t41 + t54) - g(3) * (t52 * pkin(4) - t53 * qJ(5) + t48);];
U_reg = t1;
