% Calculate inertial parameters regressor of potential energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RPPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:35
% EndTime: 2018-11-14 13:45:35
% DurationCPUTime: 0.08s
% Computational Cost: add. (176->46), mult. (198->56), div. (0->0), fcn. (177->10), ass. (0->33)
t77 = pkin(4) + pkin(6);
t83 = sin(t77) / 0.2e1;
t63 = cos(pkin(4));
t82 = t63 * qJ(2) + pkin(5);
t64 = sin(qJ(1));
t65 = cos(qJ(1));
t61 = sin(pkin(4));
t80 = qJ(2) * t61;
t81 = t65 * pkin(1) + t64 * t80;
t78 = pkin(4) - pkin(6);
t74 = sin(t78);
t79 = t83 - t74 / 0.2e1;
t76 = t65 * t80;
t75 = cos(t77);
t60 = sin(pkin(6));
t71 = cos(t78) / 0.2e1;
t68 = t71 + t75 / 0.2e1;
t42 = t64 * t60 - t65 * t68;
t62 = cos(pkin(6));
t43 = t64 * t62 + t65 * t79;
t58 = t64 * pkin(1);
t73 = t43 * pkin(2) + t42 * qJ(3) + t58;
t72 = g(1) * t64 - g(2) * t65;
t50 = t83 + t74 / 0.2e1;
t51 = t71 - t75 / 0.2e1;
t70 = t51 * pkin(2) - t50 * qJ(3) + t82;
t44 = t65 * t60 + t64 * t68;
t45 = t65 * t62 - t64 * t79;
t69 = t45 * pkin(2) + t44 * qJ(3) + t81;
t67 = g(1) * t44 + g(2) * t42 - g(3) * t50;
t66 = g(1) * t45 + g(2) * t43 + g(3) * t51;
t46 = -g(3) * t63 - t72 * t61;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t65 - g(2) * t64, t72, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t66, t67, t46, -g(1) * t81 - g(2) * (t58 - t76) - g(3) * t82, 0, 0, 0, 0, 0, 0, t46, t66, -t67, -g(1) * t69 - g(2) * (t73 - t76) - g(3) * t70, 0, 0, 0, 0, 0, 0, t46, -t67, -t66, -g(1) * (t64 * t61 * pkin(3) + t45 * qJ(4) + t69) - g(2) * (t43 * qJ(4) + (-pkin(3) - qJ(2)) * t65 * t61 + t73) - g(3) * (t63 * pkin(3) + t51 * qJ(4) + t70);];
U_reg  = t1;
