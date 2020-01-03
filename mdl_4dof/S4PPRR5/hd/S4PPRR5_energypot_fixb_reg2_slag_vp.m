% Calculate inertial parameters regressor of potential energy for
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:51
% EndTime: 2019-12-31 16:19:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (46->37), mult. (81->45), div. (0->0), fcn. (75->6), ass. (0->22)
t29 = cos(pkin(6));
t43 = g(2) * t29;
t33 = cos(qJ(3));
t42 = g(3) * t33;
t41 = g(3) * qJ(1);
t30 = sin(qJ(4));
t31 = sin(qJ(3));
t40 = t30 * t31;
t32 = cos(qJ(4));
t39 = t31 * t32;
t38 = pkin(2) + qJ(1);
t28 = sin(pkin(6));
t37 = t29 * pkin(1) + t28 * qJ(2);
t36 = t29 * pkin(4) + t37;
t25 = t28 * pkin(1);
t35 = -t29 * qJ(2) + t25;
t34 = pkin(3) * t31 - pkin(5) * t33;
t21 = g(1) * t28 - t43;
t24 = t28 * pkin(4);
t22 = g(1) * t29 + g(2) * t28;
t20 = -g(3) * t31 + t21 * t33;
t1 = [0, 0, 0, 0, 0, 0, -t22, t21, -g(3), -t41, 0, 0, 0, 0, 0, 0, -g(3), t22, -t21, -g(1) * t37 - g(2) * t35 - t41, 0, 0, 0, 0, 0, 0, -t21 * t31 - t42, -t20, -t22, -g(1) * t36 - g(2) * (t24 + t35) - g(3) * t38, 0, 0, 0, 0, 0, 0, -g(1) * (t28 * t39 + t29 * t30) - g(2) * (t28 * t30 - t29 * t39) - t32 * t42, -g(1) * (-t28 * t40 + t29 * t32) - g(2) * (t28 * t32 + t29 * t40) + t30 * t42, t20, -g(1) * (t34 * t28 + t36) - g(2) * (t24 + t25) - g(3) * (t33 * pkin(3) + t31 * pkin(5) + t38) - (-qJ(2) - t34) * t43;];
U_reg = t1;
