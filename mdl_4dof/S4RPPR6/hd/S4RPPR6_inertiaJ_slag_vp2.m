% Calculate joint inertia matrix for
% S4RPPR6
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:38
% EndTime: 2019-12-31 16:40:38
% DurationCPUTime: 0.18s
% Computational Cost: add. (109->56), mult. (222->77), div. (0->0), fcn. (165->4), ass. (0->20)
t14 = sin(pkin(6));
t15 = cos(pkin(6));
t25 = t14 ^ 2 + t15 ^ 2;
t21 = t14 * qJ(3) + pkin(1);
t6 = -t15 * pkin(2) - t21;
t24 = -0.2e1 * t6;
t23 = t25 * qJ(2) ^ 2;
t22 = -pkin(5) + qJ(2);
t16 = sin(qJ(4));
t17 = cos(qJ(4));
t4 = -t14 * t16 - t15 * t17;
t5 = t14 * t17 - t15 * t16;
t19 = -t4 * mrSges(5,1) + t5 * mrSges(5,2);
t11 = t14 * mrSges(3,2);
t8 = t22 * t15;
t7 = t22 * t14;
t3 = (pkin(2) + pkin(3)) * t15 + t21;
t2 = t16 * t7 + t17 * t8;
t1 = -t16 * t8 + t17 * t7;
t9 = [Ifges(2,3) + Ifges(5,1) * t5 ^ 2 + 0.2e1 * t3 * t19 - 0.2e1 * pkin(1) * t11 + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t24 + (Ifges(4,3) + Ifges(3,2)) * t15) * t15 + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t6 ^ 2 + t23) + m(3) * (pkin(1) ^ 2 + t23) + (0.2e1 * Ifges(5,4) * t5 + Ifges(5,2) * t4) * t4 + 0.2e1 * (-t1 * t5 + t2 * t4) * mrSges(5,3) + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * qJ(2) * t25 + (mrSges(4,3) * t24 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t15 + (Ifges(3,1) + Ifges(4,1)) * t14) * t14; -m(3) * pkin(1) - t14 * mrSges(4,3) + t11 + (-mrSges(4,1) - mrSges(3,1)) * t15 + m(4) * t6 - m(5) * t3 - t19; m(3) + m(4) + m(5); m(5) * (t17 * t1 + t16 * t2) + (m(4) * qJ(2) + mrSges(4,2)) * t14 + (t16 * t4 - t17 * t5) * mrSges(5,3); 0; m(4) + m(5) * (t16 ^ 2 + t17 ^ 2); t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t5 + Ifges(5,6) * t4; 0; t17 * mrSges(5,1) - t16 * mrSges(5,2); Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1), t9(2), t9(4), t9(7); t9(2), t9(3), t9(5), t9(8); t9(4), t9(5), t9(6), t9(9); t9(7), t9(8), t9(9), t9(10);];
Mq = res;
