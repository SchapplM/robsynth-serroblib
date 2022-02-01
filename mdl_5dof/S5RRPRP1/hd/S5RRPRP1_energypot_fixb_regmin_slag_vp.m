% Calculate minimal parameter regressor of potential energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:58
% EndTime: 2022-01-20 10:19:58
% DurationCPUTime: 0.13s
% Computational Cost: add. (70->24), mult. (51->32), div. (0->0), fcn. (45->8), ass. (0->19)
t50 = qJ(1) + qJ(2);
t46 = sin(t50);
t53 = sin(qJ(1));
t59 = t53 * pkin(1) + pkin(2) * t46;
t47 = cos(t50);
t55 = cos(qJ(1));
t58 = t55 * pkin(1) + pkin(2) * t47;
t57 = qJ(3) + pkin(6) + pkin(5);
t45 = pkin(8) + t50;
t40 = sin(t45);
t41 = cos(t45);
t56 = g(1) * t41 + g(2) * t40;
t54 = cos(qJ(4));
t52 = sin(qJ(4));
t51 = -qJ(5) - pkin(7);
t44 = pkin(4) * t54 + pkin(3);
t39 = -g(3) * t52 - t54 * t56;
t38 = -g(3) * t54 + t52 * t56;
t1 = [0, -g(1) * t55 - g(2) * t53, g(1) * t53 - g(2) * t55, 0, -g(1) * t47 - g(2) * t46, g(1) * t46 - g(2) * t47, -g(1) * t58 - g(2) * t59 - g(3) * t57, 0, 0, 0, 0, 0, t39, t38, t39, t38, -g(1) * t40 + g(2) * t41, -g(1) * (-t40 * t51 + t41 * t44 + t58) - g(2) * (t40 * t44 + t41 * t51 + t59) - g(3) * (pkin(4) * t52 + t57);];
U_reg = t1;
