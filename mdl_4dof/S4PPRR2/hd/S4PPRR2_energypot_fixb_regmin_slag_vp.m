% Calculate minimal parameter regressor of potential energy for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% 
% Output:
% U_reg [1x8]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:00
% EndTime: 2019-03-08 18:17:00
% DurationCPUTime: 0.02s
% Computational Cost: add. (19->9), mult. (12->11), div. (0->0), fcn. (8->4), ass. (0->8)
t48 = g(2) * qJ(1);
t47 = pkin(6) + qJ(3);
t46 = qJ(4) + t47;
t45 = cos(t47);
t44 = sin(t47);
t43 = cos(t46);
t42 = sin(t46);
t1 = [-t48, -g(1) * pkin(1) - g(3) * qJ(2) - t48, 0, -g(1) * t45 - g(2) * t44, g(1) * t44 - g(2) * t45, 0, -g(1) * t43 - g(2) * t42, g(1) * t42 - g(2) * t43;];
U_reg  = t1;
