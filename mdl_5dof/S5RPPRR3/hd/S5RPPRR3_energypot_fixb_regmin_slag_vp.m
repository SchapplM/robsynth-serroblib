% Calculate minimal parameter regressor of potential energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:35
% EndTime: 2022-01-23 09:14:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (57->20), mult. (45->28), div. (0->0), fcn. (39->10), ass. (0->15)
t65 = g(3) * (qJ(2) + pkin(5));
t57 = pkin(9) + qJ(4);
t58 = qJ(1) + pkin(8);
t53 = sin(t58);
t55 = cos(t58);
t64 = g(1) * t55 + g(2) * t53;
t61 = sin(qJ(1));
t62 = cos(qJ(1));
t63 = -g(1) * t62 - g(2) * t61;
t56 = qJ(5) + t57;
t54 = cos(t57);
t52 = sin(t57);
t51 = cos(t56);
t50 = sin(t56);
t1 = [0, t63, g(1) * t61 - g(2) * t62, t63 * pkin(1) - t65, -g(3) * sin(pkin(9)) - t64 * cos(pkin(9)), -g(1) * t53 + g(2) * t55, -g(1) * (t62 * pkin(1) + t55 * pkin(2) + t53 * qJ(3)) - g(2) * (t61 * pkin(1) + t53 * pkin(2) - t55 * qJ(3)) - t65, 0, 0, 0, 0, 0, -g(3) * t52 - t64 * t54, -g(3) * t54 + t64 * t52, 0, 0, 0, 0, 0, -g(3) * t50 - t64 * t51, -g(3) * t51 + t64 * t50;];
U_reg = t1;
