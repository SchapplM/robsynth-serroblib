% Calculate potential energy for
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
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S3RPP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_energypot_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_energypot_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_energypot_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPP1_energypot_fixb_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:04:50
% EndTime: 2019-03-08 18:04:50
% DurationCPUTime: 0.06s
% Computational Cost: add. (31->29), mult. (38->32), div. (0->0), fcn. (22->2), ass. (0->6)
t16 = rSges(4,3) + qJ(3);
t13 = sin(qJ(1));
t14 = cos(qJ(1));
t15 = t14 * pkin(1) + t13 * qJ(2);
t11 = t13 * pkin(1);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t14 * rSges(2,1) - t13 * rSges(2,2)) + g(2) * (t13 * rSges(2,1) + t14 * rSges(2,2)) + g(3) * (pkin(3) + rSges(2,3))) - m(3) * (g(1) * (-t14 * rSges(3,2) + t13 * rSges(3,3) + t15) + g(2) * (-t13 * rSges(3,2) + t11 + (-rSges(3,3) - qJ(2)) * t14) + g(3) * (pkin(3) + rSges(3,1))) - m(4) * (g(1) * (t13 * rSges(4,2) + t15) + g(2) * (t16 * t13 + t11) + g(3) * (pkin(2) + pkin(3) + rSges(4,1)) + (g(1) * t16 + g(2) * (-rSges(4,2) - qJ(2))) * t14);
U  = t1;
