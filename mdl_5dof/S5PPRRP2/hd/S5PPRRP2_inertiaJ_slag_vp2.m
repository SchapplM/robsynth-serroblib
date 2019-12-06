% Calculate joint inertia matrix for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:33
% EndTime: 2019-12-05 15:08:34
% DurationCPUTime: 0.16s
% Computational Cost: add. (116->57), mult. (283->68), div. (0->0), fcn. (196->6), ass. (0->24)
t21 = cos(qJ(4));
t16 = t21 ^ 2;
t19 = sin(qJ(4));
t27 = t19 ^ 2 + t16;
t36 = m(5) + m(6);
t10 = -t21 * mrSges(5,1) + t19 * mrSges(5,2);
t9 = -t21 * mrSges(6,1) - t19 * mrSges(6,3);
t35 = t10 + t9;
t34 = -m(6) * pkin(4) - mrSges(6,1);
t33 = (mrSges(6,2) + mrSges(5,3)) * t27;
t17 = sin(pkin(8));
t18 = cos(pkin(8));
t20 = sin(qJ(3));
t29 = cos(qJ(3));
t5 = t20 * t17 - t29 * t18;
t32 = t5 ^ 2;
t7 = t29 * t17 + t20 * t18;
t30 = t27 * pkin(6) * t7;
t28 = t27 * pkin(6) ^ 2;
t24 = t21 * pkin(4) + t19 * qJ(5);
t23 = (m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t21 + (-mrSges(5,1) + t34) * t19;
t8 = -pkin(3) - t24;
t4 = t7 ^ 2;
t1 = [m(2) + m(3) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t4 + t32) + (t27 * t4 + t32) * t36; 0; t27 * t36 + m(3) + m(4); (-mrSges(4,1) + t35) * t5 + m(5) * (-pkin(3) * t5 + t30) + m(6) * (t8 * t5 + t30) + (-mrSges(4,2) + t33) * t7; 0; -0.2e1 * pkin(3) * t10 + 0.2e1 * t8 * t9 + Ifges(4,3) + m(5) * (pkin(3) ^ 2 + t28) + m(6) * (t8 ^ 2 + t28) + (Ifges(5,2) + Ifges(6,3)) * t16 + ((Ifges(5,1) + Ifges(6,1)) * t19 + 0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t21) * t19 + 0.2e1 * pkin(6) * t33; t23 * t7; m(6) * t24 - t35; (qJ(5) * mrSges(6,2) + Ifges(5,6) - Ifges(6,6)) * t21 + (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t19 + t23 * pkin(6); Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); m(6) * t19 * t7; -m(6) * t21; (m(6) * pkin(6) + mrSges(6,2)) * t19; t34; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
