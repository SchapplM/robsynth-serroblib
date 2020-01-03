% Calculate inertial parameters regressor of potential energy for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:03
% EndTime: 2019-12-31 17:10:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (71->51), mult. (115->68), div. (0->0), fcn. (113->8), ass. (0->28)
t72 = g(3) * pkin(4);
t56 = sin(qJ(2));
t71 = g(3) * t56;
t55 = -pkin(6) - qJ(3);
t70 = t55 * t56;
t53 = sin(pkin(7));
t57 = sin(qJ(1));
t69 = t57 * t53;
t58 = cos(qJ(2));
t68 = t57 * t58;
t52 = pkin(7) + qJ(4);
t46 = sin(t52);
t59 = cos(qJ(1));
t67 = t59 * t46;
t47 = cos(t52);
t66 = t59 * t47;
t65 = t59 * t53;
t54 = cos(pkin(7));
t64 = t59 * t54;
t63 = t59 * pkin(1) + t57 * pkin(5);
t49 = t57 * pkin(1);
t62 = -t59 * pkin(5) + t49;
t61 = g(1) * t59 + g(2) * t57;
t60 = pkin(2) * t58 + qJ(3) * t56;
t45 = t54 * pkin(3) + pkin(2);
t44 = g(1) * t57 - g(2) * t59;
t43 = -g(3) * t58 + t61 * t56;
t1 = [0, 0, 0, 0, 0, 0, -t61, t44, -g(3), -t72, 0, 0, 0, 0, 0, 0, -t61 * t58 - t71, t43, -t44, -g(1) * t63 - g(2) * t62 - t72, 0, 0, 0, 0, 0, 0, -g(1) * (t58 * t64 + t69) - g(2) * (t54 * t68 - t65) - t54 * t71, -g(1) * (t57 * t54 - t58 * t65) - g(2) * (-t53 * t68 - t64) + t53 * t71, -t43, -g(1) * (t60 * t59 + t63) - g(2) * (t60 * t57 + t62) - g(3) * (t56 * pkin(2) - t58 * qJ(3) + pkin(4)), 0, 0, 0, 0, 0, 0, -g(1) * (t57 * t46 + t58 * t66) - g(2) * (t47 * t68 - t67) - t47 * t71, -g(1) * (t57 * t47 - t58 * t67) - g(2) * (-t46 * t68 - t66) + t46 * t71, -t43, -g(1) * (pkin(3) * t69 + t63) - g(2) * (t45 * t68 - t57 * t70 + t49) - g(3) * (t56 * t45 + t58 * t55 + pkin(4)) + (-g(1) * (t45 * t58 - t70) - g(2) * (-pkin(3) * t53 - pkin(5))) * t59;];
U_reg = t1;
