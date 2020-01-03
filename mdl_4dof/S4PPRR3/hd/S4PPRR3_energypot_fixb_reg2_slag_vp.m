% Calculate inertial parameters regressor of potential energy for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:26
% EndTime: 2019-12-31 16:17:26
% DurationCPUTime: 0.05s
% Computational Cost: add. (50->28), mult. (82->34), div. (0->0), fcn. (86->6), ass. (0->19)
t42 = g(3) * (-pkin(4) + qJ(1));
t41 = cos(qJ(3));
t40 = sin(qJ(3));
t39 = g(3) * qJ(1);
t28 = cos(pkin(6));
t37 = sin(pkin(6));
t38 = t28 * pkin(1) + t37 * qJ(2);
t36 = t28 * pkin(2) + t38;
t35 = t37 * pkin(1) - t28 * qJ(2);
t34 = t37 * pkin(2) + t35;
t16 = -t28 * t41 - t37 * t40;
t17 = t28 * t40 - t37 * t41;
t33 = g(1) * t17 - g(2) * t16;
t32 = g(1) * t16 + g(2) * t17;
t31 = cos(qJ(4));
t30 = sin(qJ(4));
t19 = -g(1) * t28 - g(2) * t37;
t18 = g(1) * t37 - g(2) * t28;
t1 = [0, 0, 0, 0, 0, 0, t19, t18, -g(3), -t39, 0, 0, 0, 0, 0, 0, t19, -g(3), -t18, -g(1) * t38 - g(2) * t35 - t39, 0, 0, 0, 0, 0, 0, t32, t33, g(3), -g(1) * t36 - g(2) * t34 - t42, 0, 0, 0, 0, 0, 0, g(3) * t30 + t32 * t31, g(3) * t31 - t32 * t30, -t33, -g(1) * (-t16 * pkin(3) + t17 * pkin(5) + t36) - g(2) * (-t17 * pkin(3) - t16 * pkin(5) + t34) - t42;];
U_reg = t1;
