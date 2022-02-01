% Calculate minimal parameter regressor of potential energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:24
% EndTime: 2022-01-20 11:17:24
% DurationCPUTime: 0.06s
% Computational Cost: add. (68->31), mult. (64->51), div. (0->0), fcn. (69->10), ass. (0->19)
t75 = g(3) * sin(pkin(9));
t63 = qJ(1) + qJ(2);
t59 = sin(t63);
t65 = cos(pkin(9));
t74 = t59 * t65;
t61 = cos(t63);
t73 = t61 * t65;
t66 = sin(qJ(4));
t72 = t65 * t66;
t68 = cos(qJ(4));
t71 = t65 * t68;
t70 = -g(1) * t61 - g(2) * t59;
t69 = cos(qJ(1));
t67 = sin(qJ(1));
t62 = qJ(4) + qJ(5);
t60 = cos(t62);
t58 = sin(t62);
t57 = g(1) * t59 - g(2) * t61;
t1 = [0, -g(1) * t69 - g(2) * t67, g(1) * t67 - g(2) * t69, 0, t70, t57, t70 * t65 - t75, -t57, -g(1) * (t69 * pkin(1) + t61 * pkin(2) + t59 * qJ(3)) - g(2) * (t67 * pkin(1) + t59 * pkin(2) - t61 * qJ(3)) - g(3) * (pkin(6) + pkin(5)), 0, 0, 0, 0, 0, -g(1) * (t59 * t66 + t61 * t71) - g(2) * (t59 * t71 - t61 * t66) - t68 * t75, -g(1) * (t59 * t68 - t61 * t72) - g(2) * (-t59 * t72 - t61 * t68) + t66 * t75, 0, 0, 0, 0, 0, -g(1) * (t59 * t58 + t60 * t73) - g(2) * (-t61 * t58 + t60 * t74) - t60 * t75, -g(1) * (-t58 * t73 + t59 * t60) - g(2) * (-t58 * t74 - t61 * t60) + t58 * t75;];
U_reg = t1;
