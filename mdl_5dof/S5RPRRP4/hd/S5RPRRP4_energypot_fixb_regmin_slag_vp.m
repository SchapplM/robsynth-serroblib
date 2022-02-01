% Calculate minimal parameter regressor of potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:32:51
% EndTime: 2022-01-23 09:32:51
% DurationCPUTime: 0.09s
% Computational Cost: add. (86->40), mult. (112->59), div. (0->0), fcn. (118->8), ass. (0->29)
t69 = qJ(3) + qJ(4);
t64 = cos(t69);
t74 = cos(qJ(3));
t60 = t74 * pkin(3) + pkin(4) * t64 + pkin(2);
t68 = -qJ(5) - pkin(7) - pkin(6);
t70 = sin(pkin(8));
t71 = cos(pkin(8));
t89 = t60 * t71 - t68 * t70;
t88 = g(3) * t70;
t63 = sin(t69);
t73 = sin(qJ(1));
t85 = t73 * t63;
t84 = t73 * t64;
t72 = sin(qJ(3));
t83 = t73 * t72;
t82 = t73 * t74;
t75 = cos(qJ(1));
t81 = t75 * t63;
t80 = t75 * t64;
t79 = t75 * t72;
t78 = t75 * t74;
t77 = t75 * pkin(1) + t73 * qJ(2);
t76 = -g(1) * t75 - g(2) * t73;
t66 = t73 * pkin(1);
t62 = g(1) * t73 - g(2) * t75;
t61 = t72 * pkin(3) + pkin(4) * t63;
t59 = -g(1) * (t71 * t80 + t85) - g(2) * (t71 * t84 - t81) - t64 * t88;
t58 = -g(1) * (-t71 * t81 + t84) - g(2) * (-t71 * t85 - t80) + t63 * t88;
t1 = [0, t76, t62, t76 * t71 - t88, -t62, -g(1) * t77 - g(2) * (-t75 * qJ(2) + t66) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(1) * (t71 * t78 + t83) - g(2) * (t71 * t82 - t79) - t74 * t88, -g(1) * (-t71 * t79 + t82) - g(2) * (-t71 * t83 - t78) + t72 * t88, 0, 0, 0, 0, 0, t59, t58, t59, t58, g(3) * t71 + t76 * t70, -g(1) * (t73 * t61 + t77) - g(2) * (t89 * t73 + t66) - g(3) * (t70 * t60 + t71 * t68 + pkin(5)) + (-g(1) * t89 - g(2) * (-qJ(2) - t61)) * t75;];
U_reg = t1;
