% Calculate potential energy for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:06
% EndTime: 2019-12-31 16:56:06
% DurationCPUTime: 0.22s
% Computational Cost: add. (61->36), mult. (96->29), div. (0->0), fcn. (75->6), ass. (0->14)
t33 = sin(qJ(4));
t36 = cos(qJ(4));
t51 = -m(5) * pkin(3) - t36 * mrSges(5,1) + t33 * mrSges(5,2) - mrSges(4,1);
t50 = -m(5) * pkin(6) + mrSges(4,2) - mrSges(5,3);
t49 = -m(4) - m(5);
t47 = -t33 * mrSges(5,1) - t36 * mrSges(5,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t34 = sin(qJ(3));
t37 = cos(qJ(3));
t46 = t51 * t34 - t50 * t37 + mrSges(2,2) - mrSges(3,3);
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t43 = t38 * pkin(1) + t35 * qJ(2);
t30 = t35 * pkin(1);
t1 = (-mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t49 * (pkin(2) + pkin(4)) + t51 * t37 + t50 * t34 + (-m(2) - m(3)) * pkin(4)) * g(3) + (-m(3) * t30 - mrSges(1,2) + t49 * (t35 * pkin(5) + t30) + ((m(3) - t49) * qJ(2) - t46) * t38 + t47 * t35) * g(2) + (-m(3) * t43 - mrSges(1,1) + t49 * (t38 * pkin(5) + t43) + t47 * t38 + t46 * t35) * g(1);
U = t1;
