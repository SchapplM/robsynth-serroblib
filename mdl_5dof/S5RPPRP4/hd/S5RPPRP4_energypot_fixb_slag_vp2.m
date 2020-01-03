% Calculate potential energy for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:00
% EndTime: 2019-12-31 17:52:00
% DurationCPUTime: 0.24s
% Computational Cost: add. (100->44), mult. (148->35), div. (0->0), fcn. (141->6), ass. (0->19)
t68 = m(5) + m(6);
t67 = mrSges(5,2) + mrSges(6,2);
t66 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t65 = -mrSges(2,1) - mrSges(3,1);
t64 = mrSges(2,2) - mrSges(3,3);
t63 = -m(4) - t68;
t49 = sin(qJ(4));
t50 = cos(qJ(4));
t62 = t68 * pkin(3) - t67 * t49 + t66 * t50 + mrSges(4,1);
t61 = m(5) * pkin(6) - m(6) * (-qJ(5) - pkin(6)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t60 = sin(qJ(1));
t51 = cos(qJ(1));
t59 = t51 * pkin(1) + t60 * qJ(2);
t58 = cos(pkin(7));
t57 = sin(pkin(7));
t55 = t60 * pkin(1) - t51 * qJ(2);
t37 = t51 * t57 - t60 * t58;
t36 = -t51 * t58 - t60 * t57;
t1 = (-mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + t67 * t50 + t66 * t49 + t63 * (-qJ(3) + pkin(5)) + (-m(2) - m(3)) * pkin(5)) * g(3) + (-m(3) * t55 - mrSges(1,2) + t65 * t60 - t64 * t51 + t63 * (t60 * pkin(2) + t55) + t62 * t37 + t61 * t36) * g(2) + (-m(3) * t59 - mrSges(1,1) + t64 * t60 + t65 * t51 + t63 * (t51 * pkin(2) + t59) - t61 * t37 + t62 * t36) * g(1);
U = t1;
