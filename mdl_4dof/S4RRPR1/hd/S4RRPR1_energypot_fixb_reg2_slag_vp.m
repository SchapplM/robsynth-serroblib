% Calculate inertial parameters regressor of potential energy for
% S4RRPR1
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
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:00
% EndTime: 2019-03-08 18:35:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (63->28), mult. (38->31), div. (0->0), fcn. (28->8), ass. (0->17)
t42 = pkin(5) + pkin(4);
t35 = qJ(1) + qJ(2);
t31 = sin(t35);
t36 = sin(qJ(1));
t41 = t36 * pkin(1) + pkin(2) * t31;
t32 = cos(t35);
t37 = cos(qJ(1));
t40 = t37 * pkin(1) + pkin(2) * t32;
t39 = qJ(3) + t42;
t30 = pkin(7) + t35;
t38 = -g(1) * t37 - g(2) * t36;
t29 = qJ(4) + t30;
t26 = cos(t30);
t25 = sin(t30);
t24 = cos(t29);
t23 = sin(t29);
t1 = [0, 0, 0, 0, 0, 0, t38, g(1) * t36 - g(2) * t37, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -g(1) * t32 - g(2) * t31, g(1) * t31 - g(2) * t32, -g(3), t38 * pkin(1) - g(3) * t42, 0, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t25, g(1) * t25 - g(2) * t26, -g(3), -g(1) * t40 - g(2) * t41 - g(3) * t39, 0, 0, 0, 0, 0, 0, -g(1) * t24 - g(2) * t23, g(1) * t23 - g(2) * t24, -g(3), -g(1) * (pkin(3) * t26 + t40) - g(2) * (pkin(3) * t25 + t41) - g(3) * (pkin(6) + t39);];
U_reg  = t1;
