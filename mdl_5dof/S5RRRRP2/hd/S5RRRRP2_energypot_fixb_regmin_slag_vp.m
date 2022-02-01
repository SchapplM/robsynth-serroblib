% Calculate minimal parameter regressor of potential energy for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:27
% EndTime: 2022-01-20 11:49:27
% DurationCPUTime: 0.04s
% Computational Cost: add. (71->24), mult. (55->29), div. (0->0), fcn. (52->8), ass. (0->17)
t55 = qJ(1) + qJ(2);
t50 = sin(t55);
t52 = cos(t55);
t60 = g(1) * t52 + g(2) * t50;
t59 = cos(qJ(1));
t58 = cos(qJ(3));
t57 = sin(qJ(1));
t56 = sin(qJ(3));
t54 = qJ(3) + qJ(4);
t53 = qJ(5) + pkin(7) + pkin(8);
t51 = cos(t54);
t49 = sin(t54);
t48 = t58 * pkin(3) + pkin(4) * t51 + pkin(2);
t47 = g(1) * t50 - g(2) * t52;
t46 = -g(3) * t49 - t60 * t51;
t45 = -g(3) * t51 + t60 * t49;
t1 = [0, -g(1) * t59 - g(2) * t57, g(1) * t57 - g(2) * t59, 0, -t60, t47, 0, 0, 0, 0, 0, -g(3) * t56 - t60 * t58, -g(3) * t58 + t60 * t56, 0, 0, 0, 0, 0, t46, t45, t46, t45, -t47, -g(1) * (t59 * pkin(1) + t52 * t48 + t53 * t50) - g(2) * (t57 * pkin(1) + t50 * t48 - t52 * t53) - g(3) * (t56 * pkin(3) + pkin(4) * t49 + pkin(5) + pkin(6));];
U_reg = t1;
