% Calculate inertial parameters regressor of potential energy for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPPRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:41
% EndTime: 2019-03-09 01:37:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (114->56), mult. (185->64), div. (0->0), fcn. (205->8), ass. (0->33)
t78 = g(3) * pkin(6);
t77 = pkin(2) + pkin(6);
t54 = qJ(4) + t77;
t76 = g(3) * t54;
t57 = sin(qJ(5));
t75 = g(3) * t57;
t74 = sin(qJ(1));
t56 = sin(qJ(6));
t59 = cos(qJ(5));
t73 = t56 * t59;
t58 = cos(qJ(6));
t72 = t58 * t59;
t60 = cos(qJ(1));
t71 = t60 * pkin(1) + t74 * qJ(2);
t70 = sin(pkin(9));
t69 = t60 * qJ(3) + t71;
t51 = t74 * pkin(1);
t68 = -t60 * qJ(2) + t51;
t67 = t74 * pkin(3) + t69;
t66 = pkin(5) * t59 + pkin(8) * t57;
t55 = cos(pkin(9));
t42 = t60 * t55 - t74 * t70;
t43 = t74 * t55 + t60 * t70;
t65 = g(1) * t43 - g(2) * t42;
t64 = g(1) * t42 + g(2) * t43;
t47 = t74 * qJ(3);
t63 = t47 + t51 + (-pkin(3) - qJ(2)) * t60;
t62 = t43 * pkin(4) - t42 * pkin(7) + t67;
t61 = -t42 * pkin(4) - t43 * pkin(7) + t63;
t45 = g(1) * t60 + g(2) * t74;
t44 = g(1) * t74 - g(2) * t60;
t39 = -g(3) * t59 + t65 * t57;
t1 = [0, 0, 0, 0, 0, 0, -t45, t44, -g(3), -t78, 0, 0, 0, 0, 0, 0, -g(3), t45, -t44, -g(1) * t71 - g(2) * t68 - t78, 0, 0, 0, 0, 0, 0, -t44, g(3), -t45, -g(1) * t69 - g(2) * (t47 + t68) - g(3) * t77, 0, 0, 0, 0, 0, 0, -t65, -t64, -g(3), -g(1) * t67 - g(2) * t63 - t76, 0, 0, 0, 0, 0, 0, -t65 * t59 - t75, t39, t64, -g(1) * t62 - g(2) * t61 - t76, 0, 0, 0, 0, 0, 0, -g(1) * (-t42 * t56 + t43 * t72) - g(2) * (-t42 * t72 - t43 * t56) - t58 * t75, -g(1) * (-t42 * t58 - t43 * t73) - g(2) * (t42 * t73 - t43 * t58) + t56 * t75, -t39, -g(1) * (t66 * t43 + t62) - g(2) * (-t66 * t42 + t61) - g(3) * (t57 * pkin(5) - t59 * pkin(8) + t54);];
U_reg  = t1;
