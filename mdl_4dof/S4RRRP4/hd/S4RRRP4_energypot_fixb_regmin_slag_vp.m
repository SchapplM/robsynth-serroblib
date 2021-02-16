% Calculate minimal parameter regressor of potential energy for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:21
% EndTime: 2021-01-15 14:30:21
% DurationCPUTime: 0.03s
% Computational Cost: add. (44->18), mult. (49->23), div. (0->0), fcn. (46->6), ass. (0->14)
t46 = sin(qJ(1));
t48 = cos(qJ(1));
t49 = g(1) * t48 + g(2) * t46;
t47 = cos(qJ(2));
t45 = sin(qJ(2));
t44 = qJ(2) + qJ(3);
t43 = -qJ(4) - pkin(6) - pkin(5);
t42 = cos(t44);
t41 = sin(t44);
t40 = g(1) * t46 - g(2) * t48;
t39 = t47 * pkin(2) + pkin(3) * t42 + pkin(1);
t38 = -g(3) * t41 - t49 * t42;
t37 = -g(3) * t42 + t49 * t41;
t1 = [0, -t49, t40, 0, 0, 0, 0, 0, -g(3) * t45 - t49 * t47, -g(3) * t47 + t49 * t45, 0, 0, 0, 0, 0, t38, t37, t38, t37, -t40, -g(1) * (t48 * t39 - t46 * t43) - g(2) * (t46 * t39 + t48 * t43) - g(3) * (t45 * pkin(2) + pkin(3) * t41 + pkin(4));];
U_reg = t1;
