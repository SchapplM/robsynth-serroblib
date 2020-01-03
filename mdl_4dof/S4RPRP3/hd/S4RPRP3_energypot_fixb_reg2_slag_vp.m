% Calculate inertial parameters regressor of potential energy for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:47
% EndTime: 2019-12-31 16:42:47
% DurationCPUTime: 0.07s
% Computational Cost: add. (63->30), mult. (59->32), div. (0->0), fcn. (49->6), ass. (0->19)
t41 = qJ(2) + pkin(4);
t48 = g(3) * t41;
t39 = qJ(1) + pkin(6);
t35 = sin(t39);
t36 = cos(t39);
t47 = g(1) * t36 + g(2) * t35;
t43 = sin(qJ(1));
t45 = cos(qJ(1));
t46 = -g(1) * t45 - g(2) * t43;
t44 = cos(qJ(3));
t42 = sin(qJ(3));
t40 = -qJ(4) - pkin(5);
t38 = t45 * pkin(1);
t37 = t43 * pkin(1);
t34 = t44 * pkin(3) + pkin(2);
t32 = g(1) * t35 - g(2) * t36;
t31 = -g(3) * t42 - t47 * t44;
t30 = -g(3) * t44 + t47 * t42;
t1 = [0, 0, 0, 0, 0, 0, t46, g(1) * t43 - g(2) * t45, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -t47, t32, -g(3), t46 * pkin(1) - t48, 0, 0, 0, 0, 0, 0, t31, t30, -t32, -g(1) * (t36 * pkin(2) + t35 * pkin(5) + t38) - g(2) * (t35 * pkin(2) - t36 * pkin(5) + t37) - t48, 0, 0, 0, 0, 0, 0, t31, t30, -t32, -g(1) * (t36 * t34 - t35 * t40 + t38) - g(2) * (t35 * t34 + t36 * t40 + t37) - g(3) * (t42 * pkin(3) + t41);];
U_reg = t1;
