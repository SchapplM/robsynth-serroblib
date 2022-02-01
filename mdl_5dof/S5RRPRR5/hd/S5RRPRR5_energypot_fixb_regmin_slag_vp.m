% Calculate minimal parameter regressor of potential energy for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:54
% EndTime: 2022-01-20 11:02:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (60->21), mult. (44->27), div. (0->0), fcn. (41->10), ass. (0->14)
t56 = pkin(9) + qJ(4);
t57 = qJ(1) + qJ(2);
t54 = sin(t57);
t55 = cos(t57);
t61 = g(1) * t55 + g(2) * t54;
t60 = cos(qJ(1));
t59 = sin(qJ(1));
t53 = qJ(5) + t56;
t52 = cos(t56);
t51 = sin(t56);
t50 = cos(t53);
t49 = sin(t53);
t48 = g(1) * t54 - g(2) * t55;
t1 = [0, -g(1) * t60 - g(2) * t59, g(1) * t59 - g(2) * t60, 0, -t61, t48, -g(3) * sin(pkin(9)) - t61 * cos(pkin(9)), -t48, -g(1) * (t60 * pkin(1) + t55 * pkin(2) + t54 * qJ(3)) - g(2) * (t59 * pkin(1) + t54 * pkin(2) - t55 * qJ(3)) - g(3) * (pkin(6) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t51 - t61 * t52, -g(3) * t52 + t61 * t51, 0, 0, 0, 0, 0, -g(3) * t49 - t61 * t50, -g(3) * t50 + t61 * t49;];
U_reg = t1;
