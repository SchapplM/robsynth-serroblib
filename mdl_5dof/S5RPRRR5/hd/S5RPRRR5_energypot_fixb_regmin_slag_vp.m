% Calculate minimal parameter regressor of potential energy for
% S5RPRRR5
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
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:49:01
% EndTime: 2022-01-20 09:49:01
% DurationCPUTime: 0.03s
% Computational Cost: add. (45->13), mult. (33->18), div. (0->0), fcn. (30->8), ass. (0->13)
t52 = qJ(1) + pkin(9) + qJ(3);
t50 = sin(t52);
t51 = cos(t52);
t61 = g(1) * t51 + g(2) * t50;
t57 = sin(qJ(1));
t59 = cos(qJ(1));
t60 = -g(1) * t59 - g(2) * t57;
t58 = cos(qJ(4));
t56 = sin(qJ(4));
t55 = qJ(4) + qJ(5);
t54 = cos(t55);
t53 = sin(t55);
t1 = [0, t60, g(1) * t57 - g(2) * t59, -g(3) * (qJ(2) + pkin(5)) + t60 * pkin(1), 0, -t61, g(1) * t50 - g(2) * t51, 0, 0, 0, 0, 0, -g(3) * t56 - t61 * t58, -g(3) * t58 + t61 * t56, 0, 0, 0, 0, 0, -g(3) * t53 - t61 * t54, -g(3) * t54 + t61 * t53;];
U_reg = t1;
