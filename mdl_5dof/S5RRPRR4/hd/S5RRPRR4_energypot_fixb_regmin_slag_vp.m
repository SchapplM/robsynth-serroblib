% Calculate minimal parameter regressor of potential energy for
% S5RRPRR4
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
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:31
% EndTime: 2022-01-20 10:48:31
% DurationCPUTime: 0.04s
% Computational Cost: add. (46->18), mult. (35->25), div. (0->0), fcn. (32->10), ass. (0->13)
t59 = qJ(1) + qJ(2);
t53 = pkin(9) + t59;
t64 = g(1) * cos(t53) + g(2) * sin(t53);
t63 = cos(qJ(1));
t62 = cos(qJ(4));
t61 = sin(qJ(1));
t60 = sin(qJ(4));
t58 = qJ(4) + qJ(5);
t57 = cos(t59);
t56 = cos(t58);
t55 = sin(t59);
t54 = sin(t58);
t1 = [0, -g(1) * t63 - g(2) * t61, g(1) * t61 - g(2) * t63, 0, -g(1) * t57 - g(2) * t55, g(1) * t55 - g(2) * t57, -g(1) * (t63 * pkin(1) + pkin(2) * t57) - g(2) * (t61 * pkin(1) + pkin(2) * t55) - g(3) * (qJ(3) + pkin(6) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t60 - t64 * t62, -g(3) * t62 + t64 * t60, 0, 0, 0, 0, 0, -g(3) * t54 - t64 * t56, -g(3) * t56 + t64 * t54;];
U_reg = t1;
