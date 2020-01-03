% Calculate inertial parameters regressor of potential energy for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:28
% EndTime: 2019-12-31 16:56:28
% DurationCPUTime: 0.06s
% Computational Cost: add. (46->37), mult. (81->43), div. (0->0), fcn. (75->6), ass. (0->24)
t52 = g(3) * pkin(4);
t51 = pkin(2) + pkin(4);
t40 = cos(qJ(1));
t50 = g(2) * t40;
t39 = cos(qJ(3));
t49 = g(3) * t39;
t35 = sin(qJ(4));
t37 = sin(qJ(1));
t48 = t37 * t35;
t38 = cos(qJ(4));
t47 = t37 * t38;
t46 = t40 * t35;
t45 = t40 * t38;
t44 = t40 * pkin(1) + t37 * qJ(2);
t43 = t40 * pkin(5) + t44;
t32 = t37 * pkin(1);
t42 = -t40 * qJ(2) + t32;
t36 = sin(qJ(3));
t41 = pkin(3) * t36 - pkin(6) * t39;
t28 = g(1) * t37 - t50;
t31 = t37 * pkin(5);
t29 = g(1) * t40 + g(2) * t37;
t27 = -g(3) * t36 + t28 * t39;
t1 = [0, 0, 0, 0, 0, 0, -t29, t28, -g(3), -t52, 0, 0, 0, 0, 0, 0, -g(3), t29, -t28, -g(1) * t44 - g(2) * t42 - t52, 0, 0, 0, 0, 0, 0, -t28 * t36 - t49, -t27, -t29, -g(1) * t43 - g(2) * (t31 + t42) - g(3) * t51, 0, 0, 0, 0, 0, 0, -g(1) * (t36 * t47 + t46) - g(2) * (-t36 * t45 + t48) - t38 * t49, -g(1) * (-t36 * t48 + t45) - g(2) * (t36 * t46 + t47) + t35 * t49, t27, -g(1) * (t41 * t37 + t43) - g(2) * (t31 + t32) - g(3) * (t39 * pkin(3) + t36 * pkin(6) + t51) - (-qJ(2) - t41) * t50;];
U_reg = t1;
