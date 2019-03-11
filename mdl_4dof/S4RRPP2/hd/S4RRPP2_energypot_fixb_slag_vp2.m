% Calculate potential energy for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:58
% EndTime: 2019-03-08 18:33:58
% DurationCPUTime: 0.08s
% Computational Cost: add. (68->33), mult. (57->24), div. (0->0), fcn. (32->4), ass. (0->11)
t37 = m(4) + m(5);
t36 = pkin(5) + pkin(4);
t34 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t32 = -m(3) * pkin(1) - mrSges(2,1);
t31 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1) - mrSges(5,1);
t29 = cos(qJ(1));
t28 = sin(qJ(1));
t27 = qJ(1) + qJ(2);
t24 = cos(t27);
t23 = sin(t27);
t1 = (-mrSges(1,3) - m(2) * pkin(4) - mrSges(2,3) - mrSges(3,3) - mrSges(4,2) - m(5) * (-qJ(4) + t36) + mrSges(5,3) + (-m(3) - m(4)) * t36) * g(3) + (-mrSges(1,2) - t29 * mrSges(2,2) + t32 * t28 + (t37 * qJ(3) - t34) * t24 + t31 * t23 - t37 * (t28 * pkin(1) + t23 * pkin(2))) * g(2) + (t28 * mrSges(2,2) + t34 * t23 + t31 * t24 + t32 * t29 - mrSges(1,1) - t37 * (t29 * pkin(1) + t24 * pkin(2) + t23 * qJ(3))) * g(1);
U  = t1;
