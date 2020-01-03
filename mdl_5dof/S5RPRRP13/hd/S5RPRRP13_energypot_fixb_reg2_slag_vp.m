% Calculate inertial parameters regressor of potential energy for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:36
% EndTime: 2019-12-31 18:59:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (79->46), mult. (145->53), div. (0->0), fcn. (144->6), ass. (0->32)
t57 = sin(qJ(3));
t60 = cos(qJ(3));
t78 = pkin(3) * t57 - pkin(7) * t60;
t77 = g(3) * pkin(5);
t76 = pkin(2) + pkin(5);
t73 = g(3) * t60;
t56 = sin(qJ(4));
t58 = sin(qJ(1));
t72 = t58 * t56;
t59 = cos(qJ(4));
t71 = t58 * t59;
t61 = cos(qJ(1));
t70 = t61 * t56;
t69 = t61 * t59;
t68 = t61 * pkin(1) + t58 * qJ(2);
t67 = t61 * pkin(6) + t68;
t66 = t60 * pkin(3) + t57 * pkin(7) + t76;
t52 = t58 * pkin(1);
t65 = -t61 * qJ(2) + t52;
t44 = g(1) * t58 - g(2) * t61;
t64 = t78 * t58 + t67;
t40 = t57 * t72 - t69;
t42 = t57 * t70 + t71;
t63 = g(1) * t40 - g(2) * t42 + t56 * t73;
t51 = t58 * pkin(6);
t62 = t51 + t52 + (-qJ(2) - t78) * t61;
t45 = g(1) * t61 + g(2) * t58;
t43 = -t57 * t69 + t72;
t41 = t57 * t71 + t70;
t39 = -g(3) * t57 + t44 * t60;
t38 = -g(1) * t41 - g(2) * t43 - t59 * t73;
t1 = [0, 0, 0, 0, 0, 0, -t45, t44, -g(3), -t77, 0, 0, 0, 0, 0, 0, -g(3), t45, -t44, -g(1) * t68 - g(2) * t65 - t77, 0, 0, 0, 0, 0, 0, -t44 * t57 - t73, -t39, -t45, -g(1) * t67 - g(2) * (t51 + t65) - g(3) * t76, 0, 0, 0, 0, 0, 0, t38, t63, t39, -g(1) * t64 - g(2) * t62 - g(3) * t66, 0, 0, 0, 0, 0, 0, t38, t39, -t63, -g(1) * (t41 * pkin(4) + t40 * qJ(5) + t64) - g(2) * (t43 * pkin(4) - t42 * qJ(5) + t62) - g(3) * ((pkin(4) * t59 + qJ(5) * t56) * t60 + t66);];
U_reg = t1;
