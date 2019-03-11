% Calculate minimal parameter regressor of potential energy for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% U_reg [1x9]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S3RPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_energypot_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:57
% EndTime: 2019-03-08 18:05:57
% DurationCPUTime: 0.02s
% Computational Cost: add. (14->11), mult. (27->19), div. (0->0), fcn. (28->4), ass. (0->9)
t40 = cos(qJ(1));
t39 = cos(qJ(3));
t38 = sin(qJ(1));
t37 = sin(qJ(3));
t36 = -g(1) * t40 - g(2) * t38;
t35 = g(1) * t38 - g(2) * t40;
t34 = -t40 * t37 + t38 * t39;
t33 = -t38 * t37 - t40 * t39;
t1 = [0, t36, t35, t36, -t35, -g(1) * (t40 * pkin(1) + t38 * qJ(2)) - g(2) * (t38 * pkin(1) - t40 * qJ(2)) - g(3) * pkin(3), 0, g(1) * t33 - g(2) * t34, -g(1) * t34 - g(2) * t33;];
U_reg  = t1;
