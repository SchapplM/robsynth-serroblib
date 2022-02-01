% Calculate minimal parameter regressor of potential energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:44
% EndTime: 2022-01-23 09:12:44
% DurationCPUTime: 0.08s
% Computational Cost: add. (90->35), mult. (96->52), div. (0->0), fcn. (95->8), ass. (0->25)
t62 = cos(qJ(4));
t50 = t62 * pkin(4) + pkin(3);
t56 = sin(pkin(8));
t57 = cos(pkin(8));
t58 = -qJ(5) - pkin(6);
t75 = t50 * t57 - t56 * t58;
t74 = g(3) * t56;
t59 = qJ(2) + pkin(5);
t73 = g(3) * t59;
t55 = qJ(1) + pkin(7);
t51 = sin(t55);
t60 = sin(qJ(4));
t71 = t51 * t60;
t69 = t57 * t60;
t68 = t57 * t62;
t61 = sin(qJ(1));
t67 = t61 * pkin(1) + t51 * pkin(2);
t52 = cos(t55);
t63 = cos(qJ(1));
t66 = t63 * pkin(1) + t52 * pkin(2) + t51 * qJ(3);
t65 = -g(1) * t52 - g(2) * t51;
t64 = -g(1) * t63 - g(2) * t61;
t46 = -g(1) * (t52 * t68 + t71) - g(2) * (t51 * t68 - t52 * t60) - t62 * t74;
t45 = -g(1) * (t51 * t62 - t52 * t69) - g(2) * (-t51 * t69 - t52 * t62) + t60 * t74;
t1 = [0, t64, g(1) * t61 - g(2) * t63, t64 * pkin(1) - t73, t65 * t57 - t74, -g(1) * t51 + g(2) * t52, -g(1) * t66 - g(2) * (-t52 * qJ(3) + t67) - t73, 0, 0, 0, 0, 0, t46, t45, t46, t45, g(3) * t57 + t65 * t56, -g(1) * (pkin(4) * t71 + t66) - g(2) * (t75 * t51 + t67) - g(3) * (t56 * t50 + t57 * t58 + t59) + (-g(1) * t75 - g(2) * (-pkin(4) * t60 - qJ(3))) * t52;];
U_reg = t1;
