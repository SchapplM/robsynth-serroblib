% Calculate inertial parameters regressor of potential energy for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:17
% EndTime: 2019-12-31 16:49:17
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->33), mult. (59->36), div. (0->0), fcn. (49->8), ass. (0->20)
t43 = qJ(2) + pkin(4);
t51 = g(3) * t43;
t41 = qJ(1) + pkin(7);
t35 = sin(t41);
t36 = cos(t41);
t50 = g(1) * t36 + g(2) * t35;
t45 = sin(qJ(1));
t47 = cos(qJ(1));
t49 = -g(1) * t47 - g(2) * t45;
t48 = -pkin(6) - pkin(5);
t46 = cos(qJ(3));
t44 = sin(qJ(3));
t42 = qJ(3) + qJ(4);
t40 = t47 * pkin(1);
t39 = t45 * pkin(1);
t38 = cos(t42);
t37 = sin(t42);
t34 = t46 * pkin(3) + pkin(2);
t32 = g(1) * t35 - g(2) * t36;
t1 = [0, 0, 0, 0, 0, 0, t49, g(1) * t45 - g(2) * t47, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -t50, t32, -g(3), t49 * pkin(1) - t51, 0, 0, 0, 0, 0, 0, -g(3) * t44 - t50 * t46, -g(3) * t46 + t50 * t44, -t32, -g(1) * (t36 * pkin(2) + t35 * pkin(5) + t40) - g(2) * (t35 * pkin(2) - t36 * pkin(5) + t39) - t51, 0, 0, 0, 0, 0, 0, -g(3) * t37 - t50 * t38, -g(3) * t38 + t50 * t37, -t32, -g(1) * (t36 * t34 - t35 * t48 + t40) - g(2) * (t35 * t34 + t36 * t48 + t39) - g(3) * (t44 * pkin(3) + t43);];
U_reg = t1;
