% Calculate minimal parameter regressor of potential energy for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:59
% EndTime: 2019-03-08 18:22:59
% DurationCPUTime: 0.02s
% Computational Cost: add. (47->18), mult. (24->20), div. (0->0), fcn. (20->6), ass. (0->9)
t49 = pkin(6) + qJ(2);
t48 = qJ(3) + t49;
t47 = cos(t49);
t46 = sin(t49);
t45 = cos(t48);
t44 = sin(t48);
t43 = -g(1) * t45 - g(2) * t44;
t42 = g(1) * t44 - g(2) * t45;
t1 = [-g(3) * qJ(1), 0, -g(1) * t47 - g(2) * t46, g(1) * t46 - g(2) * t47, 0, t43, t42, t43, -t42, -g(1) * (t45 * pkin(3) + t44 * qJ(4) + pkin(2) * t47 + cos(pkin(6)) * pkin(1)) - g(2) * (t44 * pkin(3) - t45 * qJ(4) + pkin(2) * t46 + sin(pkin(6)) * pkin(1)) - g(3) * (pkin(5) + pkin(4) + qJ(1));];
U_reg  = t1;
