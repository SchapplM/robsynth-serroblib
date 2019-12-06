% Calculate inertial parameters regressor of potential energy for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:39
% EndTime: 2019-12-05 16:19:40
% DurationCPUTime: 0.07s
% Computational Cost: add. (116->46), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t63 = pkin(5) + qJ(1);
t69 = g(3) * t63;
t65 = cos(qJ(3));
t47 = t65 * pkin(3) + pkin(2);
t62 = -pkin(6) - qJ(4);
t59 = qJ(3) + pkin(9);
t64 = sin(qJ(3));
t68 = t64 * pkin(3) + t63;
t58 = pkin(8) + qJ(2);
t48 = sin(t58);
t50 = cos(t58);
t67 = g(1) * t50 + g(2) * t48;
t60 = sin(pkin(8));
t61 = cos(pkin(8));
t66 = -g(1) * t61 - g(2) * t60;
t57 = -pkin(7) + t62;
t54 = t61 * pkin(1);
t53 = t60 * pkin(1);
t52 = qJ(5) + t59;
t51 = cos(t59);
t49 = sin(t59);
t46 = cos(t52);
t45 = sin(t52);
t43 = pkin(4) * t51 + t47;
t42 = g(1) * t48 - g(2) * t50;
t1 = [0, 0, 0, 0, 0, 0, t66, g(1) * t60 - g(2) * t61, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t67, t42, -g(3), t66 * pkin(1) - t69, 0, 0, 0, 0, 0, 0, -g(3) * t64 - t67 * t65, -g(3) * t65 + t67 * t64, -t42, -g(1) * (t50 * pkin(2) + t48 * pkin(6) + t54) - g(2) * (t48 * pkin(2) - t50 * pkin(6) + t53) - t69, 0, 0, 0, 0, 0, 0, -g(3) * t49 - t67 * t51, -g(3) * t51 + t67 * t49, -t42, -g(1) * (t50 * t47 - t48 * t62 + t54) - g(2) * (t48 * t47 + t50 * t62 + t53) - g(3) * t68, 0, 0, 0, 0, 0, 0, -g(3) * t45 - t67 * t46, -g(3) * t46 + t67 * t45, -t42, -g(1) * (t50 * t43 - t48 * t57 + t54) - g(2) * (t48 * t43 + t50 * t57 + t53) - g(3) * (pkin(4) * t49 + t68);];
U_reg = t1;
