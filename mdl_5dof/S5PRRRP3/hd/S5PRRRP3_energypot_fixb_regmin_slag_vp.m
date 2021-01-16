% Calculate minimal parameter regressor of potential energy for
% S5PRRRP3
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
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:26
% EndTime: 2021-01-15 16:23:26
% DurationCPUTime: 0.11s
% Computational Cost: add. (70->23), mult. (52->26), div. (0->0), fcn. (48->8), ass. (0->15)
t48 = pkin(8) + qJ(2);
t43 = sin(t48);
t44 = cos(t48);
t52 = g(1) * t44 + g(2) * t43;
t51 = cos(qJ(3));
t50 = sin(qJ(3));
t49 = qJ(3) + qJ(4);
t47 = -qJ(5) - pkin(7) - pkin(6);
t46 = cos(t49);
t45 = sin(t49);
t42 = t51 * pkin(3) + pkin(4) * t46 + pkin(2);
t41 = g(1) * t43 - g(2) * t44;
t40 = -g(3) * t45 - t52 * t46;
t39 = -g(3) * t46 + t52 * t45;
t1 = [-g(3) * qJ(1), 0, -t52, t41, 0, 0, 0, 0, 0, -g(3) * t50 - t52 * t51, -g(3) * t51 + t52 * t50, 0, 0, 0, 0, 0, t40, t39, t40, t39, -t41, -g(1) * (t44 * t42 - t43 * t47 + cos(pkin(8)) * pkin(1)) - g(2) * (t43 * t42 + t44 * t47 + sin(pkin(8)) * pkin(1)) - g(3) * (t50 * pkin(3) + pkin(4) * t45 + pkin(5) + qJ(1));];
U_reg = t1;
