% Calculate inertial parameters regressor of potential energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:48
% EndTime: 2019-07-18 13:28:48
% DurationCPUTime: 0.05s
% Computational Cost: add. (51->28), mult. (74->34), div. (0->0), fcn. (67->8), ass. (0->23)
t57 = cos(qJ(3));
t67 = pkin(2) * t57;
t52 = qJ(3) + qJ(4);
t49 = sin(t52);
t66 = g(2) * t49;
t54 = sin(qJ(3));
t65 = g(2) * t54;
t64 = g(3) * qJ(1);
t53 = sin(qJ(5));
t55 = sin(qJ(2));
t63 = t55 * t53;
t56 = cos(qJ(5));
t62 = t55 * t56;
t58 = cos(qJ(2));
t61 = t58 * t53;
t60 = t58 * t56;
t59 = g(1) * t58 + g(3) * t55;
t50 = cos(t52);
t48 = -g(1) * pkin(1) - t64;
t47 = g(1) * t55 - g(3) * t58;
t46 = g(2) * t50 + t59 * t49;
t45 = -g(1) * (t58 * t67 + pkin(1)) + pkin(2) * t65 - g(3) * (t55 * t67 + qJ(1));
t1 = [0, 0, 0, 0, 0, 0, -g(1), -g(2), -g(3), -t64, 0, 0, 0, 0, 0, 0, -t59, t47, g(2), t48, 0, 0, 0, 0, 0, 0, -t59 * t57 + t65, g(2) * t57 + t59 * t54, -t47, t48, 0, 0, 0, 0, 0, 0, -t59 * t50 + t66, t46, -t47, t45, 0, 0, 0, 0, 0, 0, -g(1) * (t50 * t60 + t63) + t56 * t66 - g(3) * (t50 * t62 - t61), -g(1) * (-t50 * t61 + t62) - t53 * t66 - g(3) * (-t50 * t63 - t60), -t46, t45;];
U_reg  = t1;
