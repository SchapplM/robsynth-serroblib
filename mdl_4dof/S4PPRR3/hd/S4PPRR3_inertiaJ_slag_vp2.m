% Calculate joint inertia matrix for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:22
% EndTime: 2019-12-31 16:17:22
% DurationCPUTime: 0.09s
% Computational Cost: add. (33->24), mult. (89->34), div. (0->0), fcn. (45->4), ass. (0->13)
t8 = cos(qJ(4));
t4 = t8 ^ 2;
t6 = sin(qJ(4));
t13 = t6 ^ 2 + t4;
t14 = m(5) * t13;
t12 = mrSges(5,3) * t13;
t11 = -mrSges(5,1) * t6 - mrSges(5,2) * t8;
t9 = cos(qJ(3));
t7 = sin(qJ(3));
t5 = t9 ^ 2;
t3 = t7 ^ 2;
t1 = -t8 * mrSges(5,1) + t6 * mrSges(5,2);
t2 = [m(2) + m(3) + m(4) + t14; 0; m(3) + m(4) * (t3 + t5) + m(5) * (t13 * t3 + t5); 0; (m(5) * pkin(3) + mrSges(4,1) - t1) * t9 + (pkin(5) * t14 - mrSges(4,2) + t12) * t7; Ifges(4,3) + Ifges(5,2) * t4 - 0.2e1 * pkin(3) * t1 + m(5) * (t13 * pkin(5) ^ 2 + pkin(3) ^ 2) + (Ifges(5,1) * t6 + 0.2e1 * Ifges(5,4) * t8) * t6 + 0.2e1 * pkin(5) * t12; t1; t11 * t7; Ifges(5,5) * t6 + Ifges(5,6) * t8 + t11 * pkin(5); Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
