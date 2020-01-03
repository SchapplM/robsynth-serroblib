% Calculate inertial parameters regressor of potential energy for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:48
% EndTime: 2019-12-31 16:39:48
% DurationCPUTime: 0.05s
% Computational Cost: add. (50->28), mult. (82->34), div. (0->0), fcn. (86->6), ass. (0->19)
t47 = g(3) * pkin(4);
t46 = g(3) * (-qJ(3) + pkin(4));
t45 = sin(qJ(1));
t36 = cos(qJ(1));
t44 = t36 * pkin(1) + t45 * qJ(2);
t43 = cos(pkin(6));
t42 = sin(pkin(6));
t41 = t36 * pkin(2) + t44;
t40 = t45 * pkin(1) - t36 * qJ(2);
t39 = t45 * pkin(2) + t40;
t21 = -t36 * t43 - t45 * t42;
t22 = t36 * t42 - t45 * t43;
t38 = g(1) * t22 - g(2) * t21;
t37 = g(1) * t21 + g(2) * t22;
t35 = cos(qJ(4));
t34 = sin(qJ(4));
t24 = -g(1) * t36 - g(2) * t45;
t23 = g(1) * t45 - g(2) * t36;
t1 = [0, 0, 0, 0, 0, 0, t24, t23, -g(3), -t47, 0, 0, 0, 0, 0, 0, t24, -g(3), -t23, -g(1) * t44 - g(2) * t40 - t47, 0, 0, 0, 0, 0, 0, t37, t38, g(3), -g(1) * t41 - g(2) * t39 - t46, 0, 0, 0, 0, 0, 0, g(3) * t34 + t37 * t35, g(3) * t35 - t37 * t34, -t38, -g(1) * (-t21 * pkin(3) + t22 * pkin(5) + t41) - g(2) * (-t22 * pkin(3) - t21 * pkin(5) + t39) - t46;];
U_reg = t1;
