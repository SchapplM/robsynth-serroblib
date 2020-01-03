% Calculate inertial parameters regressor of potential energy for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:13
% EndTime: 2019-12-31 17:22:13
% DurationCPUTime: 0.04s
% Computational Cost: add. (70->29), mult. (48->32), div. (0->0), fcn. (38->8), ass. (0->18)
t51 = pkin(5) + pkin(4);
t50 = g(3) * (pkin(6) + t51);
t41 = qJ(1) + qJ(2);
t35 = sin(t41);
t43 = sin(qJ(1));
t49 = t43 * pkin(1) + pkin(2) * t35;
t36 = cos(t41);
t45 = cos(qJ(1));
t48 = t45 * pkin(1) + pkin(2) * t36;
t37 = qJ(3) + t41;
t33 = sin(t37);
t34 = cos(t37);
t47 = g(1) * t34 + g(2) * t33;
t46 = -g(1) * t45 - g(2) * t43;
t44 = cos(qJ(4));
t42 = sin(qJ(4));
t30 = g(1) * t33 - g(2) * t34;
t1 = [0, 0, 0, 0, 0, 0, t46, g(1) * t43 - g(2) * t45, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -g(1) * t36 - g(2) * t35, g(1) * t35 - g(2) * t36, -g(3), t46 * pkin(1) - g(3) * t51, 0, 0, 0, 0, 0, 0, -t47, t30, -g(3), -g(1) * t48 - g(2) * t49 - t50, 0, 0, 0, 0, 0, 0, -g(3) * t42 - t47 * t44, -g(3) * t44 + t47 * t42, -t30, -g(1) * (t34 * pkin(3) + t33 * pkin(7) + t48) - g(2) * (t33 * pkin(3) - t34 * pkin(7) + t49) - t50;];
U_reg = t1;
