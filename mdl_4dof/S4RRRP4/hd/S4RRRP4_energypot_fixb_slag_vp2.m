% Calculate potential energy for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:07
% EndTime: 2019-12-31 17:15:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (77->35), mult. (86->25), div. (0->0), fcn. (61->6), ass. (0->14)
t50 = -m(4) - m(5);
t49 = mrSges(4,2) + mrSges(5,2);
t48 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t33 = qJ(2) + qJ(3);
t28 = sin(t33);
t29 = cos(t33);
t34 = sin(qJ(2));
t36 = cos(qJ(2));
t47 = -m(3) * pkin(1) - t36 * mrSges(3,1) + t34 * mrSges(3,2) - mrSges(2,1) + t50 * (t36 * pkin(2) + pkin(1)) + t48 * t29 + t49 * t28;
t38 = -pkin(6) - pkin(5);
t46 = mrSges(2,2) + m(5) * (-qJ(4) + t38) - mrSges(5,3) + m(4) * t38 - mrSges(4,3) - m(3) * pkin(5) - mrSges(3,3);
t37 = cos(qJ(1));
t35 = sin(qJ(1));
t1 = (-t34 * mrSges(3,1) - t36 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t50 * (t34 * pkin(2) + pkin(4)) - t49 * t29 + t48 * t28 + (-m(2) - m(3)) * pkin(4)) * g(3) + (t47 * t35 - t46 * t37 - mrSges(1,2)) * g(2) + (t46 * t35 + t47 * t37 - mrSges(1,1)) * g(1);
U = t1;
