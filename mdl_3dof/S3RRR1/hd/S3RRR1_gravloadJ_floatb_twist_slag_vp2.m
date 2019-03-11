% Calculate Gravitation load on the joints for
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
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [3x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S3RRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:57
% EndTime: 2019-03-08 18:07:57
% DurationCPUTime: 0.09s
% Computational Cost: add. (65->18), mult. (50->20), div. (0->0), fcn. (32->6), ass. (0->15)
t18 = m(4) * pkin(2) + mrSges(3,1);
t17 = mrSges(2,1) + (m(3) + m(4)) * pkin(1);
t8 = qJ(1) + qJ(2);
t7 = qJ(3) + t8;
t3 = sin(t7);
t4 = cos(t7);
t14 = t4 * mrSges(4,1) - t3 * mrSges(4,2);
t13 = -t3 * mrSges(4,1) - t4 * mrSges(4,2);
t5 = sin(t8);
t6 = cos(t8);
t12 = t5 * mrSges(3,2) - t18 * t6 - t14;
t11 = t6 * mrSges(3,2) + t18 * t5 - t13;
t10 = cos(qJ(1));
t9 = sin(qJ(1));
t1 = [(t9 * mrSges(2,2) - t17 * t10 + t12) * g(2) + (t10 * mrSges(2,2) + t17 * t9 + t11) * g(1), t11 * g(1) + t12 * g(2), -g(1) * t13 - g(2) * t14];
taug  = t1(:);
