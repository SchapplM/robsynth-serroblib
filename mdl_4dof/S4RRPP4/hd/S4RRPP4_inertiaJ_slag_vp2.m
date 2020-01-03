% Calculate joint inertia matrix for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:58:54
% DurationCPUTime: 0.18s
% Computational Cost: add. (98->67), mult. (194->70), div. (0->0), fcn. (86->2), ass. (0->18)
t11 = sin(qJ(2));
t12 = cos(qJ(2));
t26 = t11 ^ 2 + t12 ^ 2;
t25 = -m(4) * pkin(2) - mrSges(4,1);
t17 = t11 * qJ(3) + pkin(1);
t2 = -t12 * pkin(2) - t17;
t24 = -0.2e1 * t2;
t23 = -2 * mrSges(5,3);
t22 = pkin(2) + pkin(3);
t21 = t12 * mrSges(5,1) + t11 * mrSges(5,2);
t20 = t26 * pkin(5) ^ 2;
t19 = mrSges(4,2) - mrSges(5,3);
t18 = pkin(5) - qJ(4);
t14 = qJ(3) ^ 2;
t4 = t18 * t12;
t3 = t18 * t11;
t1 = t22 * t12 + t17;
t5 = [Ifges(2,3) + 0.2e1 * t1 * t21 + m(4) * (t2 ^ 2 + t20) + m(5) * (t1 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(3) * (pkin(1) ^ 2 + t20) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t24 + t4 * t23 + (Ifges(3,2) + Ifges(5,2) + Ifges(4,3)) * t12) * t12 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t24 + t3 * t23 + (Ifges(3,1) + Ifges(5,1) + Ifges(4,1)) * t11 + 0.2e1 * (Ifges(3,4) - Ifges(5,4) - Ifges(4,5)) * t12) * t11 + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(5) * t26; m(5) * (qJ(3) * t4 - t22 * t3) - t3 * mrSges(5,1) + t4 * mrSges(5,2) + (t19 * qJ(3) + Ifges(3,6) - Ifges(4,6) + Ifges(5,6)) * t12 + (-pkin(2) * mrSges(4,2) + mrSges(5,3) * t22 + Ifges(4,4) + Ifges(3,5) - Ifges(5,5)) * t11 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t12 + (-mrSges(3,1) + t25) * t11) * pkin(5); 0.2e1 * pkin(2) * mrSges(4,1) + 0.2e1 * t22 * mrSges(5,1) + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + 0.2e1 * (mrSges(4,3) + mrSges(5,2)) * qJ(3) + m(4) * (pkin(2) ^ 2 + t14) + m(5) * (t22 ^ 2 + t14); m(5) * t3 + (m(4) * pkin(5) + t19) * t11; -m(5) * t22 - mrSges(5,1) + t25; m(4) + m(5); m(5) * t1 + t21; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1), t5(2), t5(4), t5(7); t5(2), t5(3), t5(5), t5(8); t5(4), t5(5), t5(6), t5(9); t5(7), t5(8), t5(9), t5(10);];
Mq = res;
