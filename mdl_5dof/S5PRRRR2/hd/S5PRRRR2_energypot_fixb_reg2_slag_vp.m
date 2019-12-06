% Calculate inertial parameters regressor of potential energy for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:50
% EndTime: 2019-12-05 17:04:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (71->33), mult. (48->32), div. (0->0), fcn. (36->8), ass. (0->20)
t49 = pkin(4) + qJ(1);
t52 = g(3) * (pkin(5) + t49);
t43 = sin(qJ(2));
t38 = t43 * pkin(2);
t45 = cos(qJ(2));
t51 = t45 * pkin(2) + pkin(1);
t50 = g(3) * qJ(1);
t41 = qJ(2) + qJ(3);
t35 = sin(t41);
t48 = pkin(3) * t35 + t38;
t36 = cos(t41);
t47 = pkin(3) * t36 + t51;
t37 = qJ(4) + t41;
t33 = sin(t37);
t34 = cos(t37);
t46 = g(1) * t34 + g(2) * t33;
t44 = cos(qJ(5));
t42 = sin(qJ(5));
t30 = g(1) * t33 - g(2) * t34;
t1 = [0, 0, 0, 0, 0, 0, -g(1), -g(2), -g(3), -t50, 0, 0, 0, 0, 0, 0, -g(1) * t45 - g(2) * t43, g(1) * t43 - g(2) * t45, -g(3), -g(1) * pkin(1) - t50, 0, 0, 0, 0, 0, 0, -g(1) * t36 - g(2) * t35, g(1) * t35 - g(2) * t36, -g(3), -g(1) * t51 - g(2) * t38 - g(3) * t49, 0, 0, 0, 0, 0, 0, -t46, t30, -g(3), -g(1) * t47 - g(2) * t48 - t52, 0, 0, 0, 0, 0, 0, -g(3) * t42 - t46 * t44, -g(3) * t44 + t46 * t42, -t30, -g(1) * (t33 * pkin(6) + t47) - g(2) * (-t34 * pkin(6) + t48) - t52;];
U_reg = t1;
