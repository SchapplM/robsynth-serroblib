% Calculate minimal parameter regressor of potential energy for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% U_reg [1x9]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S3RRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energypot_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:07
% EndTime: 2019-03-08 18:08:07
% DurationCPUTime: 0.02s
% Computational Cost: add. (18->8), mult. (12->12), div. (0->0), fcn. (12->6), ass. (0->9)
t41 = qJ(1) + qJ(2);
t43 = cos(qJ(1));
t42 = sin(qJ(1));
t40 = qJ(3) + t41;
t39 = cos(t41);
t38 = sin(t41);
t37 = cos(t40);
t36 = sin(t40);
t1 = [0, -g(1) * t43 - g(2) * t42, g(1) * t42 - g(2) * t43, 0, -g(1) * t39 - g(2) * t38, g(1) * t38 - g(2) * t39, 0, -g(1) * t37 - g(2) * t36, g(1) * t36 - g(2) * t37;];
U_reg  = t1;
