% Calculate inertial parameters regressor of potential energy for
% S4RPRR1
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
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:49
% EndTime: 2019-03-08 18:31:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (63->28), mult. (38->31), div. (0->0), fcn. (28->8), ass. (0->17)
t41 = qJ(2) + pkin(4);
t34 = qJ(1) + pkin(7);
t29 = sin(t34);
t35 = sin(qJ(1));
t40 = t35 * pkin(1) + pkin(2) * t29;
t30 = cos(t34);
t36 = cos(qJ(1));
t39 = t36 * pkin(1) + pkin(2) * t30;
t38 = pkin(5) + t41;
t31 = qJ(3) + t34;
t37 = -g(1) * t36 - g(2) * t35;
t28 = qJ(4) + t31;
t27 = cos(t31);
t26 = sin(t31);
t23 = cos(t28);
t22 = sin(t28);
t1 = [0, 0, 0, 0, 0, 0, t37, g(1) * t35 - g(2) * t36, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -g(1) * t30 - g(2) * t29, g(1) * t29 - g(2) * t30, -g(3), t37 * pkin(1) - g(3) * t41, 0, 0, 0, 0, 0, 0, -g(1) * t27 - g(2) * t26, g(1) * t26 - g(2) * t27, -g(3), -g(1) * t39 - g(2) * t40 - g(3) * t38, 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t22, g(1) * t22 - g(2) * t23, -g(3), -g(1) * (pkin(3) * t27 + t39) - g(2) * (pkin(3) * t26 + t40) - g(3) * (pkin(6) + t38);];
U_reg  = t1;
