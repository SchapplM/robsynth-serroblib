% Calculate inertial parameters regressor of potential energy for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:56
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->33), mult. (59->36), div. (0->0), fcn. (49->8), ass. (0->20)
t43 = qJ(2) + pkin(4);
t49 = g(3) * t43;
t40 = qJ(1) + pkin(6);
t34 = sin(t40);
t36 = cos(t40);
t48 = g(1) * t36 + g(2) * t34;
t45 = sin(qJ(1));
t46 = cos(qJ(1));
t47 = -g(1) * t46 - g(2) * t45;
t44 = -pkin(5) - qJ(3);
t42 = cos(pkin(7));
t41 = sin(pkin(7));
t39 = pkin(7) + qJ(4);
t38 = t46 * pkin(1);
t37 = t45 * pkin(1);
t35 = cos(t39);
t33 = sin(t39);
t32 = t42 * pkin(3) + pkin(2);
t30 = g(1) * t34 - g(2) * t36;
t1 = [0, 0, 0, 0, 0, 0, t47, g(1) * t45 - g(2) * t46, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -t48, t30, -g(3), t47 * pkin(1) - t49, 0, 0, 0, 0, 0, 0, -g(3) * t41 - t48 * t42, -g(3) * t42 + t48 * t41, -t30, -g(1) * (t36 * pkin(2) + t34 * qJ(3) + t38) - g(2) * (t34 * pkin(2) - t36 * qJ(3) + t37) - t49, 0, 0, 0, 0, 0, 0, -g(3) * t33 - t48 * t35, -g(3) * t35 + t48 * t33, -t30, -g(1) * (t36 * t32 - t34 * t44 + t38) - g(2) * (t34 * t32 + t36 * t44 + t37) - g(3) * (t41 * pkin(3) + t43);];
U_reg = t1;
