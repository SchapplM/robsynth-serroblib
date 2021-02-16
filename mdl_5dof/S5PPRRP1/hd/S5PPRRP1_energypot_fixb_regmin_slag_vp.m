% Calculate minimal parameter regressor of potential energy for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:48:15
% EndTime: 2021-01-15 14:48:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (77->32), mult. (90->42), div. (0->0), fcn. (91->8), ass. (0->22)
t37 = pkin(8) + qJ(3);
t35 = sin(t37);
t52 = g(3) * t35;
t51 = g(3) * qJ(1);
t38 = sin(pkin(7));
t42 = sin(qJ(4));
t50 = t38 * t42;
t43 = cos(qJ(4));
t49 = t38 * t43;
t39 = cos(pkin(7));
t48 = t39 * t42;
t47 = t39 * t43;
t46 = pkin(4) * t42 + pkin(5) + qJ(2);
t45 = g(1) * t39 + g(2) * t38;
t34 = t43 * pkin(4) + pkin(3);
t36 = cos(t37);
t40 = -qJ(5) - pkin(6);
t44 = t34 * t36 - t35 * t40 + cos(pkin(8)) * pkin(2) + pkin(1);
t32 = -g(3) * t36 + t45 * t35;
t31 = -g(1) * (t36 * t47 + t50) - g(2) * (t36 * t49 - t48) - t43 * t52;
t30 = -g(1) * (-t36 * t48 + t49) - g(2) * (-t36 * t50 - t47) + t42 * t52;
t1 = [-t51, -g(1) * (t39 * pkin(1) + t38 * qJ(2)) - g(2) * (t38 * pkin(1) - t39 * qJ(2)) - t51, 0, -t45 * t36 - t52, t32, 0, 0, 0, 0, 0, t31, t30, t31, t30, -t32, -g(3) * (t35 * t34 + t36 * t40 + sin(pkin(8)) * pkin(2) + qJ(1)) + (-g(1) * t44 + g(2) * t46) * t39 + (-g(1) * t46 - g(2) * t44) * t38;];
U_reg = t1;
