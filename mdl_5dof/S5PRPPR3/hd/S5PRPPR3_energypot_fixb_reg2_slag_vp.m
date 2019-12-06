% Calculate inertial parameters regressor of potential energy for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:54
% EndTime: 2019-12-05 15:26:54
% DurationCPUTime: 0.14s
% Computational Cost: add. (113->52), mult. (130->62), div. (0->0), fcn. (121->8), ass. (0->32)
t45 = cos(pkin(7));
t43 = qJ(2) + pkin(8);
t39 = sin(t43);
t56 = qJ(4) * t39;
t40 = cos(t43);
t62 = t40 * t45;
t66 = pkin(3) * t62 + t45 * t56;
t65 = g(3) * t40;
t64 = g(3) * qJ(1);
t44 = sin(pkin(7));
t63 = t40 * t44;
t47 = sin(qJ(5));
t61 = t44 * t47;
t49 = cos(qJ(5));
t60 = t44 * t49;
t59 = t45 * t47;
t58 = t45 * t49;
t50 = cos(qJ(2));
t38 = t50 * pkin(2) + pkin(1);
t46 = -qJ(3) - pkin(5);
t57 = t44 * t38 + t45 * t46;
t48 = sin(qJ(2));
t55 = t48 * pkin(2) + qJ(1);
t35 = t45 * t38;
t54 = -t44 * t46 + t35;
t53 = pkin(3) * t63 + t44 * t56 + t57;
t52 = g(1) * t45 + g(2) * t44;
t51 = t39 * pkin(3) - t40 * qJ(4) + t55;
t29 = g(1) * t44 - g(2) * t45;
t28 = g(3) * t39 + t40 * t52;
t27 = t39 * t52 - t65;
t1 = [0, 0, 0, 0, 0, 0, -t52, t29, -g(3), -t64, 0, 0, 0, 0, 0, 0, -g(3) * t48 - t50 * t52, -g(3) * t50 + t48 * t52, -t29, -g(1) * (t45 * pkin(1) + t44 * pkin(5)) - g(2) * (t44 * pkin(1) - t45 * pkin(5)) - t64, 0, 0, 0, 0, 0, 0, -t28, t27, -t29, -g(1) * t54 - g(2) * t57 - g(3) * t55, 0, 0, 0, 0, 0, 0, -t29, t28, -t27, -g(1) * (t54 + t66) - g(2) * t53 - g(3) * t51, 0, 0, 0, 0, 0, 0, -g(1) * (t39 * t59 + t60) - g(2) * (t39 * t61 - t58) + t47 * t65, -g(1) * (t39 * t58 - t61) - g(2) * (t39 * t60 + t59) + t49 * t65, -t28, -g(1) * (pkin(6) * t62 + t35 + (pkin(4) - t46) * t44 + t66) - g(2) * (-t45 * pkin(4) + pkin(6) * t63 + t53) - g(3) * (t39 * pkin(6) + t51);];
U_reg = t1;
