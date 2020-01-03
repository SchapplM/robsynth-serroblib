% Calculate inertial parameters regressor of potential energy for
% S5RPRRP9
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
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:10
% EndTime: 2019-12-31 18:49:10
% DurationCPUTime: 0.07s
% Computational Cost: add. (119->43), mult. (106->50), div. (0->0), fcn. (93->8), ass. (0->25)
t79 = g(3) * pkin(5);
t68 = sin(pkin(8));
t78 = t68 * pkin(2) + pkin(5);
t69 = cos(pkin(8));
t59 = t69 * pkin(2) + pkin(1);
t70 = -pkin(6) - qJ(2);
t67 = pkin(8) + qJ(3);
t61 = cos(t67);
t53 = pkin(3) * t61 + t59;
t66 = -pkin(7) + t70;
t71 = sin(qJ(1));
t72 = cos(qJ(1));
t77 = t71 * t53 + t72 * t66;
t60 = sin(t67);
t76 = pkin(3) * t60 + t78;
t75 = t72 * t53 - t71 * t66;
t74 = g(1) * t72 + g(2) * t71;
t62 = qJ(4) + t67;
t57 = sin(t62);
t58 = cos(t62);
t73 = pkin(4) * t58 + qJ(5) * t57;
t54 = g(1) * t71 - g(2) * t72;
t50 = -g(3) * t57 - t74 * t58;
t49 = -g(3) * t58 + t74 * t57;
t1 = [0, 0, 0, 0, 0, 0, -t74, t54, -g(3), -t79, 0, 0, 0, 0, 0, 0, -g(3) * t68 - t74 * t69, -g(3) * t69 + t74 * t68, -t54, -g(1) * (t72 * pkin(1) + t71 * qJ(2)) - g(2) * (t71 * pkin(1) - t72 * qJ(2)) - t79, 0, 0, 0, 0, 0, 0, -g(3) * t60 - t74 * t61, -g(3) * t61 + t74 * t60, -t54, -g(1) * (t72 * t59 - t71 * t70) - g(2) * (t71 * t59 + t72 * t70) - g(3) * t78, 0, 0, 0, 0, 0, 0, t50, t49, -t54, -g(1) * t75 - g(2) * t77 - g(3) * t76, 0, 0, 0, 0, 0, 0, t50, -t54, -t49, -g(1) * (t73 * t72 + t75) - g(2) * (t73 * t71 + t77) - g(3) * (t57 * pkin(4) - t58 * qJ(5) + t76);];
U_reg = t1;
