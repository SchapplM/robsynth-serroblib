% Calculate inertial parameters regressor of potential energy for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:09
% EndTime: 2019-12-31 17:13:09
% DurationCPUTime: 0.07s
% Computational Cost: add. (63->30), mult. (59->32), div. (0->0), fcn. (49->6), ass. (0->19)
t48 = pkin(5) + pkin(4);
t51 = g(3) * t48;
t42 = qJ(1) + qJ(2);
t38 = sin(t42);
t39 = cos(t42);
t50 = g(1) * t39 + g(2) * t38;
t45 = sin(qJ(1));
t47 = cos(qJ(1));
t49 = -g(1) * t47 - g(2) * t45;
t46 = cos(qJ(3));
t44 = sin(qJ(3));
t43 = -qJ(4) - pkin(6);
t41 = t47 * pkin(1);
t40 = t45 * pkin(1);
t37 = t46 * pkin(3) + pkin(2);
t35 = g(1) * t38 - g(2) * t39;
t34 = -g(3) * t44 - t50 * t46;
t33 = -g(3) * t46 + t50 * t44;
t1 = [0, 0, 0, 0, 0, 0, t49, g(1) * t45 - g(2) * t47, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -t50, t35, -g(3), t49 * pkin(1) - t51, 0, 0, 0, 0, 0, 0, t34, t33, -t35, -g(1) * (t39 * pkin(2) + t38 * pkin(6) + t41) - g(2) * (t38 * pkin(2) - t39 * pkin(6) + t40) - t51, 0, 0, 0, 0, 0, 0, t34, t33, -t35, -g(1) * (t39 * t37 - t38 * t43 + t41) - g(2) * (t38 * t37 + t39 * t43 + t40) - g(3) * (t44 * pkin(3) + t48);];
U_reg = t1;
