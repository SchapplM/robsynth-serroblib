% Calculate inertial parameters regressor of potential energy for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:29
% EndTime: 2019-12-31 17:33:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (76->34), mult. (121->39), div. (0->0), fcn. (132->6), ass. (0->22)
t36 = -pkin(5) + qJ(1);
t50 = g(3) * t36;
t49 = cos(qJ(3));
t48 = sin(qJ(3));
t47 = g(3) * qJ(1);
t35 = cos(pkin(7));
t45 = sin(pkin(7));
t46 = t35 * pkin(1) + t45 * qJ(2);
t44 = t35 * pkin(2) + t46;
t43 = t45 * pkin(1) - t35 * qJ(2);
t42 = t45 * pkin(2) + t43;
t23 = -t35 * t49 - t45 * t48;
t24 = t35 * t48 - t45 * t49;
t41 = g(1) * t24 - g(2) * t23;
t20 = g(1) * t23 + g(2) * t24;
t40 = -t23 * pkin(3) + t24 * qJ(4) + t44;
t39 = -t24 * pkin(3) - t23 * qJ(4) + t42;
t38 = cos(qJ(5));
t37 = sin(qJ(5));
t26 = -g(1) * t35 - g(2) * t45;
t25 = g(1) * t45 - g(2) * t35;
t1 = [0, 0, 0, 0, 0, 0, t26, t25, -g(3), -t47, 0, 0, 0, 0, 0, 0, t26, -g(3), -t25, -g(1) * t46 - g(2) * t43 - t47, 0, 0, 0, 0, 0, 0, t20, t41, g(3), -g(1) * t44 - g(2) * t42 - t50, 0, 0, 0, 0, 0, 0, g(3), -t20, -t41, -g(1) * t40 - g(2) * t39 - t50, 0, 0, 0, 0, 0, 0, g(3) * t38 - t41 * t37, -g(3) * t37 - t41 * t38, t20, -g(1) * (-t23 * pkin(6) + t40) - g(2) * (-t24 * pkin(6) + t39) - g(3) * (-pkin(4) + t36);];
U_reg = t1;
