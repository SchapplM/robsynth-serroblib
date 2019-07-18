% Calculate inertial parameters regressor of potential energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:25
% EndTime: 2019-07-18 17:22:25
% DurationCPUTime: 0.07s
% Computational Cost: add. (71->36), mult. (100->50), div. (0->0), fcn. (95->8), ass. (0->29)
t64 = cos(qJ(2));
t79 = pkin(1) * t64;
t58 = qJ(2) + qJ(4);
t55 = sin(t58);
t78 = pkin(4) * t55;
t77 = g(3) * t55;
t61 = sin(qJ(2));
t76 = g(3) * t61;
t66 = pkin(2) + pkin(1);
t75 = t61 * t66;
t60 = sin(qJ(5));
t62 = sin(qJ(1));
t74 = t62 * t60;
t63 = cos(qJ(5));
t73 = t62 * t63;
t72 = t64 * t66;
t65 = cos(qJ(1));
t71 = t65 * t60;
t70 = t65 * t63;
t59 = -pkin(3) - qJ(3);
t69 = t65 * t59 + t62 * t72;
t68 = -t62 * t59 + t65 * t72;
t67 = g(1) * t65 + g(2) * t62;
t56 = cos(t58);
t51 = g(1) * t62 - g(2) * t65;
t50 = -t67 * t64 - t76;
t49 = -g(3) * t64 + t67 * t61;
t48 = -g(3) * t56 + t67 * t55;
t1 = [0, 0, 0, 0, 0, 0, -t67, t51, -g(3), 0, 0, 0, 0, 0, 0, 0, t50, t49, -t51, 0, 0, 0, 0, 0, 0, 0, t50, t49, -t51, -g(1) * (t62 * qJ(3) + t65 * t79) - g(2) * (-t65 * qJ(3) + t62 * t79) - pkin(1) * t76, 0, 0, 0, 0, 0, 0, -t67 * t56 - t77, t48, -t51, -g(1) * t68 - g(2) * t69 - g(3) * t75, 0, 0, 0, 0, 0, 0, -g(1) * (t56 * t70 + t74) - g(2) * (t56 * t73 - t71) - t63 * t77, -g(1) * (-t56 * t71 + t73) - g(2) * (-t56 * t74 - t70) + t60 * t77, -t48, -g(1) * (t65 * t78 + t68) - g(2) * (t62 * t78 + t69) - g(3) * (-t56 * pkin(4) + t75);];
U_reg  = t1;
