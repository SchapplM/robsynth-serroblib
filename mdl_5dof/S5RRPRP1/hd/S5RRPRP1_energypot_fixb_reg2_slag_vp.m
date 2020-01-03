% Calculate inertial parameters regressor of potential energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:11
% EndTime: 2020-01-03 11:59:11
% DurationCPUTime: 0.06s
% Computational Cost: add. (112->38), mult. (74->41), div. (0->0), fcn. (61->8), ass. (0->23)
t64 = pkin(6) + pkin(5);
t52 = qJ(3) + t64;
t63 = g(1) * t52;
t53 = qJ(1) + qJ(2);
t49 = sin(t53);
t56 = sin(qJ(1));
t62 = t56 * pkin(1) + pkin(2) * t49;
t50 = cos(t53);
t58 = cos(qJ(1));
t61 = -t58 * pkin(1) - pkin(2) * t50;
t48 = pkin(8) + t53;
t44 = sin(t48);
t45 = cos(t48);
t60 = g(2) * t44 - g(3) * t45;
t59 = -g(2) * t56 + g(3) * t58;
t57 = cos(qJ(4));
t55 = sin(qJ(4));
t54 = -qJ(5) - pkin(7);
t47 = t57 * pkin(4) + pkin(3);
t43 = g(2) * t45 + g(3) * t44;
t42 = -g(1) * t57 + t60 * t55;
t41 = -g(1) * t55 - t60 * t57;
t1 = [0, 0, 0, 0, 0, 0, t59, -g(2) * t58 - g(3) * t56, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, -g(2) * t49 + g(3) * t50, -g(2) * t50 - g(3) * t49, -g(1), t59 * pkin(1) - g(1) * t64, 0, 0, 0, 0, 0, 0, -t60, -t43, -g(1), -g(2) * t62 - g(3) * t61 - t63, 0, 0, 0, 0, 0, 0, t41, t42, t43, -t63 - g(2) * (t44 * pkin(3) - t45 * pkin(7) + t62) - g(3) * (-t45 * pkin(3) - t44 * pkin(7) + t61), 0, 0, 0, 0, 0, 0, t41, t42, t43, -g(1) * (t55 * pkin(4) + t52) - g(2) * (t44 * t47 + t45 * t54 + t62) - g(3) * (t44 * t54 - t45 * t47 + t61);];
U_reg = t1;
