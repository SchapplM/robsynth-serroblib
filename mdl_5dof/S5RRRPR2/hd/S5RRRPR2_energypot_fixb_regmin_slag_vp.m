% Calculate minimal parameter regressor of potential energy for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:52
% EndTime: 2022-01-20 11:30:52
% DurationCPUTime: 0.05s
% Computational Cost: add. (49->21), mult. (31->27), div. (0->0), fcn. (28->10), ass. (0->13)
t43 = qJ(1) + qJ(2);
t42 = qJ(3) + t43;
t37 = pkin(9) + t42;
t48 = g(1) * cos(t37) + g(2) * sin(t37);
t47 = cos(qJ(1));
t46 = cos(qJ(5));
t45 = sin(qJ(1));
t44 = sin(qJ(5));
t41 = cos(t43);
t40 = sin(t43);
t39 = cos(t42);
t38 = sin(t42);
t1 = [0, -g(1) * t47 - g(2) * t45, g(1) * t45 - g(2) * t47, 0, -g(1) * t41 - g(2) * t40, g(1) * t40 - g(2) * t41, 0, -g(1) * t39 - g(2) * t38, g(1) * t38 - g(2) * t39, -g(1) * (t47 * pkin(1) + pkin(2) * t41 + pkin(3) * t39) - g(2) * (t45 * pkin(1) + pkin(2) * t40 + pkin(3) * t38) - g(3) * (qJ(4) + pkin(7) + pkin(6) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t44 - t48 * t46, -g(3) * t46 + t48 * t44;];
U_reg = t1;
