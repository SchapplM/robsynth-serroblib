% Calculate inertial parameters regressor of potential energy for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% 
% Output:
% U_reg [1x(2*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S2RR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_energypot_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_energypot_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:26
% EndTime: 2020-06-19 09:14:26
% DurationCPUTime: 0.05s
% Computational Cost: add. (14->10), mult. (14->11), div. (0->0), fcn. (10->4), ass. (0->7)
t10 = sin(qJ(1));
t11 = cos(qJ(1));
t12 = -g(1) * t11 - g(2) * t10;
t9 = qJ(1) + qJ(2);
t8 = cos(t9);
t7 = sin(t9);
t1 = [0, 0, 0, 0, 0, 0, t12, g(1) * t10 - g(2) * t11, -g(3), -g(3) * pkin(2), 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t7, g(1) * t7 - g(2) * t8, -g(3), -g(3) * (pkin(3) + pkin(2)) + t12 * pkin(1);];
U_reg = t1;
