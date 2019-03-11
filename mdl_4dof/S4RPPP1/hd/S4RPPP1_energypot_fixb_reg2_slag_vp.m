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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
% StartTime: 2019-03-08 18:26:25
% EndTime: 2019-03-08 18:26:25
% DurationCPUTime: 0.07s
% Computational Cost: add. (77->40), mult. (165->53), div. (0->0), fcn. (177->6), ass. (0->31)
t58 = sin(qJ(1));
t76 = g(1) * t58;
t59 = cos(qJ(1));
t75 = g(2) * t59;
t57 = cos(pkin(4));
t74 = t57 * qJ(2) + pkin(5);
t54 = sin(pkin(6));
t55 = sin(pkin(4));
t73 = t54 * t55;
t72 = t58 * t54;
t56 = cos(pkin(6));
t71 = t58 * t56;
t70 = t59 * t54;
t69 = t59 * t56;
t67 = qJ(2) * t55;
t68 = t59 * pkin(1) + t58 * t67;
t66 = qJ(3) * t56;
t65 = pkin(2) * t73 + t74;
t64 = t59 * t67;
t43 = -t57 * t69 + t72;
t44 = t57 * t70 + t71;
t52 = t58 * pkin(1);
t63 = t44 * pkin(2) + t43 * qJ(3) + t52;
t62 = -t75 + t76;
t45 = t57 * t71 + t70;
t46 = -t57 * t72 + t69;
t61 = t46 * pkin(2) + t45 * qJ(3) + t68;
t38 = g(3) * t55 * t56 - g(1) * t45 - g(2) * t43;
t60 = g(1) * t46 + g(2) * t44 + g(3) * t73;
t40 = -g(3) * t57 - t62 * t55;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t59 - g(2) * t58, t62, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t60, -t38, t40, -g(1) * t68 - g(2) * (t52 - t64) - g(3) * t74, 0, 0, 0, 0, 0, 0, t40, t60, t38, -g(1) * t61 - g(2) * (t63 - t64) - g(3) * (-t55 * t66 + t65) 0, 0, 0, 0, 0, 0, t40, t38, -t60, -g(1) * (t46 * qJ(4) + t61) - g(2) * (t44 * qJ(4) + t63) - g(3) * (t57 * pkin(3) + t65) + (-pkin(3) * t76 - g(3) * (qJ(4) * t54 - t66) - (-pkin(3) - qJ(2)) * t75) * t55;];
U_reg  = t1;
