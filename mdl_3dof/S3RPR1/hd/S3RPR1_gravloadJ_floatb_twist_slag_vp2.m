% Calculate Gravitation load on the joints for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S3RPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:51
% EndTime: 2019-03-08 18:05:51
% DurationCPUTime: 0.10s
% Computational Cost: add. (35->21), mult. (62->26), div. (0->0), fcn. (54->4), ass. (0->12)
t19 = m(3) + m(4);
t17 = cos(qJ(3));
t16 = sin(qJ(3));
t15 = mrSges(2,1) + mrSges(3,1);
t14 = mrSges(2,2) - mrSges(3,3);
t10 = cos(qJ(1));
t9 = sin(qJ(1));
t1 = -t10 * t17 - t9 * t16;
t2 = t10 * t16 - t9 * t17;
t12 = -t2 * mrSges(4,1) + t1 * mrSges(4,2);
t11 = t1 * mrSges(4,1) + t2 * mrSges(4,2);
t3 = [(t14 * t9 + (-m(4) * pkin(2) - t15) * t10 + t11 - t19 * (t10 * pkin(1) + t9 * qJ(2))) * g(2) + ((m(3) * pkin(1) - m(4) * (-pkin(1) - pkin(2)) + t15) * t9 + t12 + (-t19 * qJ(2) + t14) * t10) * g(1), t19 * (-g(1) * t9 + g(2) * t10) -g(1) * t12 - g(2) * t11];
taug  = t3(:);
