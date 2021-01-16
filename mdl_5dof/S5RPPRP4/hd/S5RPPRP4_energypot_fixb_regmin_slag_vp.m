% Calculate minimal parameter regressor of potential energy for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:13
% EndTime: 2021-01-15 17:13:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (58->30), mult. (88->41), div. (0->0), fcn. (89->6), ass. (0->22)
t44 = sin(pkin(7));
t45 = cos(pkin(7));
t46 = qJ(5) + pkin(6);
t49 = cos(qJ(4));
t54 = pkin(4) * t49 + pkin(3);
t57 = t54 * t44 - t46 * t45 + qJ(2);
t51 = pkin(1) + pkin(2);
t56 = -qJ(3) + pkin(5);
t50 = cos(qJ(1));
t55 = t50 * qJ(2);
t48 = sin(qJ(1));
t38 = t50 * t44 - t48 * t45;
t39 = t48 * t44 + t50 * t45;
t52 = g(1) * t39 - g(2) * t38;
t47 = sin(qJ(4));
t43 = t48 * qJ(2);
t41 = -g(1) * t50 - g(2) * t48;
t40 = g(1) * t48 - g(2) * t50;
t37 = t46 * t44 + t54 * t45 + t51;
t36 = g(3) * t47 - t52 * t49;
t35 = g(3) * t49 + t52 * t47;
t1 = [0, t41, t40, t41, -t40, -g(1) * (t50 * pkin(1) + t43) - g(2) * (t48 * pkin(1) - t55) - g(3) * pkin(5), -g(1) * (t51 * t50 + t43) - g(2) * (t51 * t48 - t55) - g(3) * t56, 0, 0, 0, 0, 0, t36, t35, t36, t35, -g(1) * t38 - g(2) * t39, -g(1) * (t37 * t50 + t57 * t48) - g(2) * (t37 * t48 - t57 * t50) - g(3) * (-t47 * pkin(4) + t56);];
U_reg = t1;
