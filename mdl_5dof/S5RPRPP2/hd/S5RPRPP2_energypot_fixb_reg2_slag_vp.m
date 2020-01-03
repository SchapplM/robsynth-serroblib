% Calculate inertial parameters regressor of potential energy for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:11
% EndTime: 2019-12-31 18:11:11
% DurationCPUTime: 0.07s
% Computational Cost: add. (109->40), mult. (104->43), div. (0->0), fcn. (91->6), ass. (0->24)
t54 = qJ(1) + pkin(7);
t48 = sin(t54);
t56 = sin(qJ(3));
t66 = qJ(4) * t56;
t58 = cos(qJ(3));
t69 = t48 * t58;
t71 = pkin(3) * t69 + t48 * t66;
t55 = qJ(2) + pkin(5);
t70 = g(3) * t55;
t49 = cos(t54);
t68 = t49 * t58;
t57 = sin(qJ(1));
t67 = t57 * pkin(1) + t48 * pkin(2);
t59 = cos(qJ(1));
t65 = t59 * pkin(1) + t49 * pkin(2) + t48 * pkin(6);
t64 = -t49 * pkin(6) + t67;
t63 = pkin(3) * t68 + t49 * t66 + t65;
t62 = g(1) * t49 + g(2) * t48;
t61 = -g(1) * t59 - g(2) * t57;
t60 = t56 * pkin(3) - t58 * qJ(4) + t55;
t39 = g(1) * t48 - g(2) * t49;
t38 = -g(3) * t56 - t62 * t58;
t37 = -g(3) * t58 + t62 * t56;
t1 = [0, 0, 0, 0, 0, 0, t61, g(1) * t57 - g(2) * t59, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t62, t39, -g(3), t61 * pkin(1) - t70, 0, 0, 0, 0, 0, 0, t38, t37, -t39, -g(1) * t65 - g(2) * t64 - t70, 0, 0, 0, 0, 0, 0, t38, -t39, -t37, -g(1) * t63 - g(2) * (t64 + t71) - g(3) * t60, 0, 0, 0, 0, 0, 0, t38, -t37, t39, -g(1) * (pkin(4) * t68 - t48 * qJ(5) + t63) - g(2) * (pkin(4) * t69 + (-pkin(6) + qJ(5)) * t49 + t67 + t71) - g(3) * (t56 * pkin(4) + t60);];
U_reg = t1;
