% Calculate inertial parameters regressor of potential energy for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:19
% EndTime: 2019-12-31 17:52:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (81->38), mult. (130->43), div. (0->0), fcn. (141->6), ass. (0->24)
t67 = g(3) * pkin(5);
t53 = -qJ(3) + pkin(5);
t66 = g(3) * t53;
t65 = sin(qJ(1));
t56 = cos(qJ(1));
t64 = t56 * pkin(1) + t65 * qJ(2);
t63 = cos(pkin(7));
t62 = sin(pkin(7));
t61 = t56 * pkin(2) + t64;
t60 = t65 * pkin(1) - t56 * qJ(2);
t59 = t65 * pkin(2) + t60;
t39 = -t56 * t63 - t65 * t62;
t40 = t56 * t62 - t65 * t63;
t58 = g(1) * t40 - g(2) * t39;
t57 = g(1) * t39 + g(2) * t40;
t55 = cos(qJ(4));
t54 = sin(qJ(4));
t52 = -qJ(5) - pkin(6);
t46 = t55 * pkin(4) + pkin(3);
t42 = -g(1) * t56 - g(2) * t65;
t41 = g(1) * t65 - g(2) * t56;
t37 = g(3) * t54 + t57 * t55;
t36 = g(3) * t55 - t57 * t54;
t1 = [0, 0, 0, 0, 0, 0, t42, t41, -g(3), -t67, 0, 0, 0, 0, 0, 0, t42, -g(3), -t41, -g(1) * t64 - g(2) * t60 - t67, 0, 0, 0, 0, 0, 0, t57, t58, g(3), -g(1) * t61 - g(2) * t59 - t66, 0, 0, 0, 0, 0, 0, t37, t36, -t58, -g(1) * (-t39 * pkin(3) + t40 * pkin(6) + t61) - g(2) * (-t40 * pkin(3) - t39 * pkin(6) + t59) - t66, 0, 0, 0, 0, 0, 0, t37, t36, -t58, -g(1) * (-t39 * t46 - t40 * t52 + t61) - g(2) * (t39 * t52 - t40 * t46 + t59) - g(3) * (-t54 * pkin(4) + t53);];
U_reg = t1;
