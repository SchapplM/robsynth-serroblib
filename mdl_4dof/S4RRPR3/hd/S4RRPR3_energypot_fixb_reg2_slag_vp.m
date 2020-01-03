% Calculate inertial parameters regressor of potential energy for
% S4RRPR3
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:33
% EndTime: 2019-12-31 17:01:33
% DurationCPUTime: 0.04s
% Computational Cost: add. (70->29), mult. (48->32), div. (0->0), fcn. (38->8), ass. (0->18)
t47 = pkin(5) + pkin(4);
t46 = g(3) * (qJ(3) + t47);
t37 = qJ(1) + qJ(2);
t32 = sin(t37);
t39 = sin(qJ(1));
t45 = t39 * pkin(1) + pkin(2) * t32;
t33 = cos(t37);
t41 = cos(qJ(1));
t44 = t41 * pkin(1) + pkin(2) * t33;
t31 = pkin(7) + t37;
t27 = sin(t31);
t28 = cos(t31);
t43 = g(1) * t28 + g(2) * t27;
t42 = -g(1) * t41 - g(2) * t39;
t40 = cos(qJ(4));
t38 = sin(qJ(4));
t26 = g(1) * t27 - g(2) * t28;
t1 = [0, 0, 0, 0, 0, 0, t42, g(1) * t39 - g(2) * t41, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -g(1) * t33 - g(2) * t32, g(1) * t32 - g(2) * t33, -g(3), t42 * pkin(1) - g(3) * t47, 0, 0, 0, 0, 0, 0, -t43, t26, -g(3), -g(1) * t44 - g(2) * t45 - t46, 0, 0, 0, 0, 0, 0, -g(3) * t38 - t43 * t40, -g(3) * t40 + t43 * t38, -t26, -g(1) * (t28 * pkin(3) + t27 * pkin(6) + t44) - g(2) * (t27 * pkin(3) - t28 * pkin(6) + t45) - t46;];
U_reg = t1;
