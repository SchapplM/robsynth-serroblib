% Calculate minimal parameter regressor of potential energy for
% S5RRPPR2
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:50
% EndTime: 2022-01-20 10:05:50
% DurationCPUTime: 0.06s
% Computational Cost: add. (71->29), mult. (53->43), div. (0->0), fcn. (51->10), ass. (0->18)
t66 = g(3) * (qJ(3) + pkin(6) + pkin(5));
t65 = g(3) * sin(pkin(9));
t56 = cos(pkin(9));
t57 = sin(qJ(5));
t64 = t56 * t57;
t59 = cos(qJ(5));
t63 = t56 * t59;
t54 = qJ(1) + qJ(2);
t49 = sin(t54);
t58 = sin(qJ(1));
t62 = t58 * pkin(1) + pkin(2) * t49;
t50 = cos(t54);
t60 = cos(qJ(1));
t61 = t60 * pkin(1) + pkin(2) * t50;
t48 = pkin(8) + t54;
t45 = cos(t48);
t44 = sin(t48);
t1 = [0, -g(1) * t60 - g(2) * t58, g(1) * t58 - g(2) * t60, 0, -g(1) * t50 - g(2) * t49, g(1) * t49 - g(2) * t50, -g(1) * t61 - g(2) * t62 - t66, -t65 + (-g(1) * t45 - g(2) * t44) * t56, -g(1) * t44 + g(2) * t45, -g(1) * (t45 * pkin(3) + t44 * qJ(4) + t61) - g(2) * (t44 * pkin(3) - t45 * qJ(4) + t62) - t66, 0, 0, 0, 0, 0, -g(1) * (t44 * t57 + t45 * t63) - g(2) * (t44 * t63 - t45 * t57) - t59 * t65, -g(1) * (t44 * t59 - t45 * t64) - g(2) * (-t44 * t64 - t45 * t59) + t57 * t65;];
U_reg = t1;
