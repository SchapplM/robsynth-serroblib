% Calculate minimal parameter regressor of potential energy for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:18
% EndTime: 2021-01-15 12:27:18
% DurationCPUTime: 0.04s
% Computational Cost: add. (53->24), mult. (62->30), div. (0->0), fcn. (56->6), ass. (0->14)
t43 = qJ(3) + qJ(4);
t40 = sin(t43);
t44 = sin(qJ(3));
t48 = t44 * pkin(3) + pkin(4) * t40 + qJ(2);
t45 = sin(qJ(1));
t47 = cos(qJ(1));
t37 = g(1) * t45 - g(2) * t47;
t46 = cos(qJ(3));
t42 = qJ(5) + pkin(1) + pkin(6) + pkin(7);
t41 = cos(t43);
t38 = g(1) * t47 + g(2) * t45;
t36 = g(3) * t40 - t37 * t41;
t35 = -g(3) * t41 - t37 * t40;
t1 = [0, -t38, t37, t38, -t37, -g(1) * (t47 * pkin(1) + t45 * qJ(2)) - g(2) * (t45 * pkin(1) - t47 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t46 - t37 * t44, g(3) * t44 - t37 * t46, 0, 0, 0, 0, 0, t35, t36, t35, t36, -t38, -g(1) * (t42 * t47 + t48 * t45) - g(2) * (t42 * t45 - t48 * t47) - g(3) * (t46 * pkin(3) + pkin(4) * t41 + pkin(2) + pkin(5));];
U_reg = t1;
