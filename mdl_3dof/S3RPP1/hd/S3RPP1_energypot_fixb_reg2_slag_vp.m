% Calculate inertial parameters regressor of potential energy for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% 
% Output:
% U_reg [1x(3*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S3RPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_energypot_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_energypot_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:01
% EndTime: 2019-03-08 18:05:01
% DurationCPUTime: 0.03s
% Computational Cost: add. (21->19), mult. (29->16), div. (0->0), fcn. (22->2), ass. (0->8)
t16 = g(3) * pkin(3);
t12 = sin(qJ(1));
t13 = cos(qJ(1));
t15 = t13 * pkin(1) + t12 * qJ(2);
t14 = t12 * pkin(1) - t13 * qJ(2);
t7 = g(1) * t13 + g(2) * t12;
t6 = g(1) * t12 - g(2) * t13;
t1 = [0, 0, 0, 0, 0, 0, -t7, t6, -g(3), -t16, 0, 0, 0, 0, 0, 0, -g(3), t7, -t6, -g(1) * t15 - g(2) * t14 - t16, 0, 0, 0, 0, 0, 0, -g(3), -t6, -t7, -g(1) * (t13 * qJ(3) + t15) - g(2) * (t12 * qJ(3) + t14) - g(3) * (pkin(2) + pkin(3));];
U_reg  = t1;
