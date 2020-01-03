% Calculate inertial parameters regressor of potential energy for
% S4RRPR4
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:32
% EndTime: 2019-12-31 17:02:33
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->33), mult. (59->36), div. (0->0), fcn. (49->8), ass. (0->20)
t50 = pkin(5) + pkin(4);
t53 = g(3) * t50;
t44 = qJ(1) + qJ(2);
t39 = sin(t44);
t40 = cos(t44);
t52 = g(1) * t40 + g(2) * t39;
t48 = sin(qJ(1));
t49 = cos(qJ(1));
t51 = -g(1) * t49 - g(2) * t48;
t47 = -pkin(6) - qJ(3);
t46 = cos(pkin(7));
t45 = sin(pkin(7));
t43 = pkin(7) + qJ(4);
t42 = t49 * pkin(1);
t41 = t48 * pkin(1);
t38 = cos(t43);
t37 = sin(t43);
t35 = t46 * pkin(3) + pkin(2);
t34 = g(1) * t39 - g(2) * t40;
t1 = [0, 0, 0, 0, 0, 0, t51, g(1) * t48 - g(2) * t49, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -t52, t34, -g(3), t51 * pkin(1) - t53, 0, 0, 0, 0, 0, 0, -g(3) * t45 - t52 * t46, -g(3) * t46 + t52 * t45, -t34, -g(1) * (t40 * pkin(2) + t39 * qJ(3) + t42) - g(2) * (t39 * pkin(2) - t40 * qJ(3) + t41) - t53, 0, 0, 0, 0, 0, 0, -g(3) * t37 - t52 * t38, -g(3) * t38 + t52 * t37, -t34, -g(1) * (t40 * t35 - t39 * t47 + t42) - g(2) * (t39 * t35 + t40 * t47 + t41) - g(3) * (t45 * pkin(3) + t50);];
U_reg = t1;
