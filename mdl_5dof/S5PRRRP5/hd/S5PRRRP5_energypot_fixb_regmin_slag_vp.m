% Calculate minimal parameter regressor of potential energy for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:41
% EndTime: 2021-01-15 16:33:41
% DurationCPUTime: 0.07s
% Computational Cost: add. (82->35), mult. (105->54), div. (0->0), fcn. (113->8), ass. (0->23)
t54 = sin(qJ(2));
t64 = g(3) * t54;
t50 = qJ(3) + qJ(4);
t47 = sin(t50);
t53 = sin(qJ(3));
t63 = t53 * pkin(3) + pkin(4) * t47 + pkin(5);
t51 = sin(pkin(8));
t56 = cos(qJ(2));
t62 = t51 * t56;
t52 = cos(pkin(8));
t61 = t52 * t56;
t60 = t53 * t56;
t55 = cos(qJ(3));
t59 = t55 * t56;
t58 = g(1) * t52 + g(2) * t51;
t48 = cos(t50);
t45 = t55 * pkin(3) + pkin(4) * t48 + pkin(2);
t49 = -qJ(5) - pkin(7) - pkin(6);
t57 = t45 * t56 - t49 * t54 + pkin(1);
t44 = -g(3) * t56 + t58 * t54;
t43 = -g(1) * (t51 * t47 + t48 * t61) - g(2) * (-t52 * t47 + t48 * t62) - t48 * t64;
t42 = -g(1) * (-t47 * t61 + t51 * t48) - g(2) * (-t47 * t62 - t52 * t48) + t47 * t64;
t1 = [-g(3) * qJ(1), 0, -t58 * t56 - t64, t44, 0, 0, 0, 0, 0, -g(1) * (t51 * t53 + t52 * t59) - g(2) * (t51 * t59 - t52 * t53) - t55 * t64, -g(1) * (t51 * t55 - t52 * t60) - g(2) * (-t51 * t60 - t52 * t55) + t53 * t64, 0, 0, 0, 0, 0, t43, t42, t43, t42, -t44, -g(3) * (t54 * t45 + t56 * t49 + qJ(1)) + (-g(1) * t57 + g(2) * t63) * t52 + (-g(1) * t63 - g(2) * t57) * t51;];
U_reg = t1;
