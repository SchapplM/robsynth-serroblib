% Calculate inertial parameters regressor of potential energy for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% U_reg [1x(2*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S2RR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_energypot_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_energypot_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:02
% EndTime: 2019-03-08 18:00:02
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->5), mult. (20->9), div. (0->0), fcn. (18->4), ass. (0->7)
t13 = sin(qJ(1));
t15 = cos(qJ(1));
t11 = g(1) * t15 - g(3) * t13;
t16 = g(1) * t13 + g(3) * t15;
t14 = cos(qJ(2));
t12 = sin(qJ(2));
t1 = [0, 0, 0, 0, 0, 0, t16, t11, -g(2), 0, 0, 0, 0, 0, 0, 0, g(2) * t12 + t16 * t14, g(2) * t14 - t16 * t12, t11, t11 * pkin(1);];
U_reg  = t1;
