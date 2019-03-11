% Calculate minimal parameter regressor of potential energy for
% S4PPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta3]';
% 
% Output:
% U_reg [1x6]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:11:16
% EndTime: 2019-03-08 18:11:16
% DurationCPUTime: 0.02s
% Computational Cost: add. (13->10), mult. (11->9), div. (0->0), fcn. (4->2), ass. (0->6)
t31 = g(1) * qJ(2);
t30 = g(2) * qJ(1);
t29 = pkin(5) + qJ(4);
t28 = cos(t29);
t27 = sin(t29);
t1 = [-t30, g(3) * pkin(1) - t30 - t31, -t31 - g(2) * (pkin(2) + qJ(1)) - g(3) * (-qJ(3) - pkin(1)) 0, -g(1) * t27 - g(2) * t28, -g(1) * t28 + g(2) * t27;];
U_reg  = t1;
