% Calculate joint inertia matrix for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR4_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:49
% DurationCPUTime: 0.10s
% Computational Cost: add. (59->31), mult. (109->42), div. (0->0), fcn. (47->4), ass. (0->15)
t10 = sin(pkin(6));
t5 = t10 * pkin(1) + qJ(3);
t18 = t5 ^ 2;
t13 = cos(qJ(4));
t9 = t13 ^ 2;
t12 = sin(qJ(4));
t16 = t12 ^ 2 + t9;
t1 = m(5) * t16;
t17 = m(4) + t1;
t15 = t16 * mrSges(5,3);
t11 = cos(pkin(6));
t7 = -t11 * pkin(1) - pkin(2);
t4 = -pkin(5) + t7;
t2 = -t12 * mrSges(5,1) - t13 * mrSges(5,2);
t3 = [Ifges(5,1) * t9 + 0.2e1 * t7 * mrSges(4,2) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (-0.2e1 * Ifges(5,4) * t13 + Ifges(5,2) * t12) * t12 + m(5) * (t16 * t4 ^ 2 + t18) + m(4) * (t7 ^ 2 + t18) + m(3) * (t10 ^ 2 + t11 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (-t2 + mrSges(4,3)) * t5 + 0.2e1 * (t11 * mrSges(3,1) - t10 * mrSges(3,2)) * pkin(1) - 0.2e1 * t4 * t15; 0; m(3) + t17; m(4) * t7 + t4 * t1 + mrSges(4,2) - t15; 0; t17; (mrSges(5,1) * t4 + Ifges(5,5)) * t13 + (-mrSges(5,2) * t4 - Ifges(5,6)) * t12; t2; t13 * mrSges(5,1) - t12 * mrSges(5,2); Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
