% Calculate inertial parameters regressor of potential energy for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR16_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:31
% EndTime: 2019-12-31 18:39:31
% DurationCPUTime: 0.08s
% Computational Cost: add. (69->50), mult. (119->58), div. (0->0), fcn. (110->6), ass. (0->29)
t77 = g(3) * pkin(5);
t76 = pkin(2) + pkin(5);
t61 = cos(qJ(1));
t75 = g(2) * t61;
t57 = sin(qJ(3));
t74 = g(3) * t57;
t58 = sin(qJ(1));
t73 = t57 * t58;
t60 = cos(qJ(3));
t72 = t58 * t60;
t56 = sin(qJ(5));
t71 = t61 * t56;
t59 = cos(qJ(5));
t70 = t61 * t59;
t69 = t61 * pkin(1) + t58 * qJ(2);
t68 = qJ(4) * t60;
t51 = t58 * pkin(6);
t52 = t58 * pkin(1);
t67 = t61 * t68 + t51 + t52;
t66 = t61 * pkin(6) + t69;
t65 = t60 * pkin(3) + t57 * qJ(4) + t76;
t64 = -pkin(3) * t57 - qJ(2);
t63 = -t61 * qJ(2) + t52;
t44 = g(1) * t58 - t75;
t62 = pkin(3) * t73 - t58 * t68 + t66;
t45 = g(1) * t61 + g(2) * t58;
t43 = t44 * t60 - t74;
t42 = g(1) * t73 + g(3) * t60 - t57 * t75;
t1 = [0, 0, 0, 0, 0, 0, -t45, t44, -g(3), -t77, 0, 0, 0, 0, 0, 0, -g(3), t45, -t44, -g(1) * t69 - g(2) * t63 - t77, 0, 0, 0, 0, 0, 0, -t42, -t43, -t45, -g(1) * t66 - g(2) * (t51 + t63) - g(3) * t76, 0, 0, 0, 0, 0, 0, -t45, t42, t43, -g(1) * t62 - g(2) * (t64 * t61 + t67) - g(3) * t65, 0, 0, 0, 0, 0, 0, -g(1) * (-t56 * t72 + t70) - g(2) * (t58 * t59 + t60 * t71) - t56 * t74, -g(1) * (-t59 * t72 - t71) - g(2) * (-t58 * t56 + t60 * t70) - t59 * t74, -t42, -g(1) * (pkin(7) * t73 + t62) - g(2) * (t58 * pkin(4) + t67) - g(3) * (t60 * pkin(7) + t65) + (-g(1) * pkin(4) - g(2) * (-pkin(7) * t57 + t64)) * t61;];
U_reg = t1;
