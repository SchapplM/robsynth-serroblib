% Calculate minimal parameter regressor of potential energy for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:41
% EndTime: 2022-01-20 09:51:41
% DurationCPUTime: 0.05s
% Computational Cost: add. (65->24), mult. (43->31), div. (0->0), fcn. (37->10), ass. (0->16)
t60 = g(3) * (qJ(3) + pkin(6) + pkin(5));
t53 = qJ(1) + qJ(2);
t47 = sin(t53);
t55 = sin(qJ(1));
t59 = t55 * pkin(1) + pkin(2) * t47;
t48 = cos(t53);
t56 = cos(qJ(1));
t58 = t56 * pkin(1) + pkin(2) * t48;
t46 = pkin(8) + t53;
t40 = sin(t46);
t41 = cos(t46);
t57 = g(1) * t41 + g(2) * t40;
t52 = pkin(9) + qJ(5);
t45 = cos(t52);
t44 = sin(t52);
t1 = [0, -g(1) * t56 - g(2) * t55, g(1) * t55 - g(2) * t56, 0, -g(1) * t48 - g(2) * t47, g(1) * t47 - g(2) * t48, -g(1) * t58 - g(2) * t59 - t60, -g(3) * sin(pkin(9)) - t57 * cos(pkin(9)), -g(1) * t40 + g(2) * t41, -g(1) * (t41 * pkin(3) + t40 * qJ(4) + t58) - g(2) * (t40 * pkin(3) - t41 * qJ(4) + t59) - t60, 0, 0, 0, 0, 0, -g(3) * t44 - t57 * t45, -g(3) * t45 + t57 * t44;];
U_reg = t1;
