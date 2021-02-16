% Calculate minimal parameter regressor of potential energy for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:51
% EndTime: 2021-01-15 15:04:51
% DurationCPUTime: 0.07s
% Computational Cost: add. (92->34), mult. (92->48), div. (0->0), fcn. (93->8), ass. (0->22)
t53 = cos(qJ(4));
t43 = t53 * pkin(4) + pkin(3);
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t51 = -qJ(5) - pkin(6);
t64 = t43 * t50 - t49 * t51;
t63 = g(3) * t49;
t48 = pkin(7) + qJ(2);
t44 = sin(t48);
t52 = sin(qJ(4));
t61 = t44 * t52;
t59 = t50 * t52;
t58 = t50 * t53;
t57 = pkin(5) + qJ(1);
t56 = t44 * pkin(2) + sin(pkin(7)) * pkin(1);
t45 = cos(t48);
t55 = t45 * pkin(2) + t44 * qJ(3) + cos(pkin(7)) * pkin(1);
t54 = -g(1) * t45 - g(2) * t44;
t39 = g(1) * t44 - g(2) * t45;
t38 = -g(1) * (t45 * t58 + t61) - g(2) * (t44 * t58 - t45 * t52) - t53 * t63;
t37 = -g(1) * (t44 * t53 - t45 * t59) - g(2) * (-t44 * t59 - t45 * t53) + t52 * t63;
t1 = [-g(3) * qJ(1), 0, t54, t39, t54 * t50 - t63, -t39, -g(1) * t55 - g(2) * (-t45 * qJ(3) + t56) - g(3) * t57, 0, 0, 0, 0, 0, t38, t37, t38, t37, g(3) * t50 + t54 * t49, -g(1) * (pkin(4) * t61 + t55) - g(2) * (t64 * t44 + t56) - g(3) * (t49 * t43 + t50 * t51 + t57) + (-g(1) * t64 - g(2) * (-pkin(4) * t52 - qJ(3))) * t45;];
U_reg = t1;
