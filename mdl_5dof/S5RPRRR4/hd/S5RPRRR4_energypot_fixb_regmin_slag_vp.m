% Calculate minimal parameter regressor of potential energy for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:45
% EndTime: 2022-01-23 09:34:45
% DurationCPUTime: 0.03s
% Computational Cost: add. (45->13), mult. (27->18), div. (0->0), fcn. (24->8), ass. (0->13)
t40 = qJ(1) + pkin(9) + qJ(3);
t39 = qJ(4) + t40;
t35 = sin(t39);
t36 = cos(t39);
t46 = g(1) * t36 + g(2) * t35;
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t45 = -g(1) * t44 - g(2) * t42;
t43 = cos(qJ(5));
t41 = sin(qJ(5));
t38 = cos(t40);
t37 = sin(t40);
t1 = [0, t45, g(1) * t42 - g(2) * t44, -g(3) * (qJ(2) + pkin(5)) + t45 * pkin(1), 0, -g(1) * t38 - g(2) * t37, g(1) * t37 - g(2) * t38, 0, -t46, g(1) * t35 - g(2) * t36, 0, 0, 0, 0, 0, -g(3) * t41 - t46 * t43, -g(3) * t43 + t46 * t41;];
U_reg = t1;
