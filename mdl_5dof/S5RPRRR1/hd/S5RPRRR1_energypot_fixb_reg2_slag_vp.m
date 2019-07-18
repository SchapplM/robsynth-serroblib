% Calculate inertial parameters regressor of potential energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:01
% EndTime: 2019-07-18 13:26:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (43->34), mult. (107->45), div. (0->0), fcn. (115->8), ass. (0->23)
t50 = sin(qJ(3));
t64 = g(3) * t50;
t51 = sin(qJ(1));
t63 = t50 * t51;
t53 = cos(qJ(4));
t62 = t50 * t53;
t55 = cos(qJ(1));
t61 = t50 * t55;
t54 = cos(qJ(3));
t60 = t51 * t54;
t49 = sin(qJ(4));
t59 = t55 * t49;
t58 = t55 * t53;
t57 = g(1) * t55 + g(2) * t51;
t45 = g(1) * t51 - g(2) * t55;
t56 = g(1) * (-t51 * t53 + t54 * t59) + g(2) * (t49 * t60 + t58) + t49 * t64;
t52 = cos(qJ(5));
t48 = sin(qJ(5));
t44 = t45 * qJ(2);
t43 = t51 * t49 + t54 * t58;
t41 = t53 * t60 - t59;
t39 = -g(3) * t54 + t57 * t50;
t1 = [0, 0, 0, 0, 0, 0, -t57, t45, -g(3), 0, 0, 0, 0, 0, 0, 0, -t57, -g(3), -t45, -t44, 0, 0, 0, 0, 0, 0, -t57 * t54 - t64, t39, -t45, -t44, 0, 0, 0, 0, 0, 0, -g(1) * t43 - g(2) * t41 - g(3) * t62, t56, -t39, -t44, 0, 0, 0, 0, 0, 0, -g(1) * (t43 * t52 + t48 * t61) - g(2) * (t41 * t52 + t48 * t63) - g(3) * (-t54 * t48 + t52 * t62), -g(1) * (-t43 * t48 + t52 * t61) - g(2) * (-t41 * t48 + t52 * t63) - g(3) * (-t48 * t62 - t54 * t52), -t56, -t44;];
U_reg  = t1;
