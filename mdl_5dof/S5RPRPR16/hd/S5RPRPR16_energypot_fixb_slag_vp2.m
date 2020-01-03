% Calculate potential energy for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR16_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:40
% EndTime: 2019-12-31 18:38:40
% DurationCPUTime: 0.29s
% Computational Cost: add. (88->55), mult. (137->45), div. (0->0), fcn. (110->6), ass. (0->19)
t69 = -mrSges(4,1) + mrSges(5,2);
t68 = m(5) + m(6);
t46 = sin(qJ(5));
t49 = cos(qJ(5));
t67 = -t46 * mrSges(6,1) - t49 * mrSges(6,2) - mrSges(5,3);
t66 = m(6) * pkin(7) + mrSges(6,3);
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t65 = -t50 * mrSges(4,2) + t69 * t47 + mrSges(2,2) - mrSges(3,3);
t64 = -m(6) * pkin(4) - t49 * mrSges(6,1) + t46 * mrSges(6,2) - mrSges(2,1) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t63 = pkin(2) + pkin(5);
t62 = pkin(3) * t47;
t48 = sin(qJ(1));
t42 = t48 * pkin(1);
t60 = t48 * pkin(6) + t42;
t51 = cos(qJ(1));
t59 = t51 * pkin(1) + t48 * qJ(2);
t57 = t51 * pkin(6) + t59;
t1 = (-m(4) * t63 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) - t68 * (t50 * pkin(3) + t47 * qJ(4) + t63) + (-m(2) - m(3)) * pkin(5) + (-t66 + t69) * t50 + (mrSges(4,2) + t67) * t47) * g(3) + (-m(3) * t42 - m(4) * t60 - mrSges(1,2) - t68 * (t51 * t50 * qJ(4) + t60) + (m(5) * t62 - (m(6) * (-pkin(3) - pkin(7)) - mrSges(6,3)) * t47 + t67 * t50 + (m(3) + m(4) + t68) * qJ(2) - t65) * t51 + t64 * t48) * g(2) + (-m(3) * t59 - m(4) * t57 - mrSges(1,1) - t68 * (t48 * t62 + t57) + t64 * t51 + (-t66 * t47 + (t68 * qJ(4) - t67) * t50 + t65) * t48) * g(1);
U = t1;
