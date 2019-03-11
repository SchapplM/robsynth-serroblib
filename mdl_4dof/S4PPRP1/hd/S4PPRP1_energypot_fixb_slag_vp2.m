% Calculate potential energy for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:13
% EndTime: 2019-03-08 18:12:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (60->36), mult. (85->29), div. (0->0), fcn. (72->4), ass. (0->14)
t50 = -m(4) - m(5);
t49 = cos(qJ(3));
t48 = sin(qJ(3));
t47 = -mrSges(2,1) - mrSges(3,1);
t46 = mrSges(2,2) - mrSges(3,3);
t38 = sin(pkin(5));
t39 = cos(pkin(5));
t44 = t39 * pkin(1) + t38 * qJ(2);
t42 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t41 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3);
t35 = t38 * pkin(1);
t29 = -t38 * t49 + t39 * t48;
t28 = -t38 * t48 - t39 * t49;
t1 = (-mrSges(3,2) + mrSges(5,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + t50 * (-pkin(4) + qJ(1)) + (-m(2) - m(3)) * qJ(1)) * g(3) + (-mrSges(1,2) - m(3) * t35 + ((m(3) - t50) * qJ(2) - t46) * t39 + t47 * t38 + t42 * t29 + t41 * t28 + t50 * (t38 * pkin(2) + t35)) * g(2) + (-m(3) * t44 + t42 * t28 - t41 * t29 + t46 * t38 + t47 * t39 - mrSges(1,1) + t50 * (t39 * pkin(2) + t44)) * g(1);
U  = t1;
