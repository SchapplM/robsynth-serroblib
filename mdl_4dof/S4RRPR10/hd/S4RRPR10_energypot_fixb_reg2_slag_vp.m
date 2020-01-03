% Calculate inertial parameters regressor of potential energy for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:55
% EndTime: 2019-12-31 17:11:55
% DurationCPUTime: 0.07s
% Computational Cost: add. (53->41), mult. (102->49), div. (0->0), fcn. (96->6), ass. (0->26)
t51 = sin(qJ(1));
t50 = sin(qJ(2));
t59 = qJ(3) * t50;
t53 = cos(qJ(2));
t64 = t51 * t53;
t69 = pkin(2) * t64 + t51 * t59;
t68 = g(3) * pkin(4);
t67 = g(3) * t53;
t49 = sin(qJ(4));
t66 = t51 * t49;
t52 = cos(qJ(4));
t65 = t51 * t52;
t54 = cos(qJ(1));
t63 = t53 * t54;
t62 = t54 * t49;
t61 = t54 * t52;
t60 = t54 * pkin(1) + t51 * pkin(5);
t46 = t51 * pkin(1);
t58 = -t54 * pkin(5) + t46;
t57 = pkin(2) * t63 + t54 * t59 + t60;
t56 = t50 * pkin(2) - t53 * qJ(3) + pkin(4);
t55 = g(1) * t54 + g(2) * t51;
t39 = g(1) * t51 - g(2) * t54;
t38 = g(3) * t50 + t55 * t53;
t37 = t55 * t50 - t67;
t1 = [0, 0, 0, 0, 0, 0, -t55, t39, -g(3), -t68, 0, 0, 0, 0, 0, 0, -t38, t37, -t39, -g(1) * t60 - g(2) * t58 - t68, 0, 0, 0, 0, 0, 0, -t39, t38, -t37, -g(1) * t57 - g(2) * (t58 + t69) - g(3) * t56, 0, 0, 0, 0, 0, 0, -g(1) * (t50 * t62 + t65) - g(2) * (t50 * t66 - t61) + t49 * t67, -g(1) * (t50 * t61 - t66) - g(2) * (t50 * t65 + t62) + t52 * t67, -t38, -g(1) * (t51 * pkin(3) + pkin(6) * t63 + t57) - g(2) * (pkin(6) * t64 + t46 + (-pkin(3) - pkin(5)) * t54 + t69) - g(3) * (t50 * pkin(6) + t56);];
U_reg = t1;
