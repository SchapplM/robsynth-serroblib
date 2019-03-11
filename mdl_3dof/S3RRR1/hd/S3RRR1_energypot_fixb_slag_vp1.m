% Calculate potential energy for
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
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S3RRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energypot_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_energypot_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRR1_energypot_fixb_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:47
% EndTime: 2019-03-08 18:07:47
% DurationCPUTime: 0.05s
% Computational Cost: add. (43->30), mult. (34->32), div. (0->0), fcn. (18->6), ass. (0->12)
t25 = pkin(3) + pkin(4);
t22 = qJ(1) + qJ(2);
t24 = cos(qJ(1));
t23 = sin(qJ(1));
t21 = t24 * pkin(1);
t20 = t23 * pkin(1);
t19 = qJ(3) + t22;
t18 = cos(t22);
t17 = sin(t22);
t16 = cos(t19);
t15 = sin(t19);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t24 * rSges(2,1) - t23 * rSges(2,2)) + g(2) * (t23 * rSges(2,1) + t24 * rSges(2,2)) + g(3) * (pkin(3) + rSges(2,3))) - m(3) * (g(1) * (t18 * rSges(3,1) - t17 * rSges(3,2) + t21) + g(2) * (t17 * rSges(3,1) + t18 * rSges(3,2) + t20) + g(3) * (rSges(3,3) + t25)) - m(4) * (g(1) * (t16 * rSges(4,1) - t15 * rSges(4,2) + pkin(2) * t18 + t21) + g(2) * (t15 * rSges(4,1) + t16 * rSges(4,2) + pkin(2) * t17 + t20) + g(3) * (pkin(5) + rSges(4,3) + t25));
U  = t1;
