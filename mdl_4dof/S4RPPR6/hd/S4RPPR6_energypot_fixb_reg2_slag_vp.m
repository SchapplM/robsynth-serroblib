% Calculate inertial parameters regressor of potential energy for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:43
% EndTime: 2019-12-31 16:40:43
% DurationCPUTime: 0.07s
% Computational Cost: add. (54->35), mult. (105->43), div. (0->0), fcn. (101->6), ass. (0->23)
t47 = sin(qJ(1));
t44 = sin(pkin(6));
t56 = qJ(3) * t44;
t45 = cos(pkin(6));
t59 = t45 * t47;
t63 = pkin(2) * t59 + t47 * t56;
t49 = cos(qJ(1));
t52 = g(1) * t49 + g(2) * t47;
t62 = g(3) * pkin(4);
t58 = t45 * t49;
t57 = t49 * pkin(1) + t47 * qJ(2);
t41 = t47 * pkin(1);
t55 = -t49 * qJ(2) + t41;
t54 = pkin(2) * t58 + t49 * t56 + t57;
t53 = t44 * pkin(2) - t45 * qJ(3) + pkin(4);
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t51 = t44 * t48 - t45 * t46;
t50 = t44 * t46 + t45 * t48;
t34 = g(1) * t47 - g(2) * t49;
t33 = -g(3) * t44 - t52 * t45;
t32 = -g(3) * t45 + t52 * t44;
t1 = [0, 0, 0, 0, 0, 0, -t52, t34, -g(3), -t62, 0, 0, 0, 0, 0, 0, t33, t32, -t34, -g(1) * t57 - g(2) * t55 - t62, 0, 0, 0, 0, 0, 0, t33, -t34, -t32, -g(1) * t54 - g(2) * (t55 + t63) - g(3) * t53, 0, 0, 0, 0, 0, 0, -g(3) * t51 - t52 * t50, g(3) * t50 - t52 * t51, t34, -g(1) * (pkin(3) * t58 - t47 * pkin(5) + t54) - g(2) * (pkin(3) * t59 + t41 + (pkin(5) - qJ(2)) * t49 + t63) - g(3) * (t44 * pkin(3) + t53);];
U_reg = t1;
