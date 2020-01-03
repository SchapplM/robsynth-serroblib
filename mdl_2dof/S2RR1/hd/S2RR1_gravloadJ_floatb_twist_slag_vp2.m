% Calculate Gravitation load on the joints for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S2RR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_gravloadJ_floatb_twist_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [3x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:06
% EndTime: 2020-01-03 11:19:07
% DurationCPUTime: 0.08s
% Computational Cost: add. (16->10), mult. (35->15), div. (0->0), fcn. (26->4), ass. (0->8)
t8 = m(3) * pkin(1) + mrSges(2,2) + mrSges(3,3);
t1 = sin(qJ(2));
t3 = cos(qJ(2));
t7 = -mrSges(3,1) * t3 + mrSges(3,2) * t1;
t5 = mrSges(2,1) - t7;
t4 = cos(qJ(1));
t2 = sin(qJ(1));
t6 = [(t5 * t2 + t8 * t4) * g(3) + (t8 * t2 - t5 * t4) * g(1), -g(2) * t7 + (t2 * g(1) + t4 * g(3)) * (mrSges(3,1) * t1 + mrSges(3,2) * t3)];
taug = t6(:);
