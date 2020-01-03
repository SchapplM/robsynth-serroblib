% Calculate joint inertia matrix for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:27
% EndTime: 2019-12-31 17:34:27
% DurationCPUTime: 0.16s
% Computational Cost: add. (90->59), mult. (212->70), div. (0->0), fcn. (117->4), ass. (0->22)
t27 = -2 * mrSges(6,3);
t26 = m(5) + m(6);
t25 = m(6) * pkin(4);
t24 = -qJ(5) - pkin(6);
t14 = sin(qJ(4));
t16 = cos(qJ(4));
t23 = t14 ^ 2 + t16 ^ 2;
t21 = mrSges(6,1) + t25;
t8 = t14 * mrSges(6,2);
t3 = -t16 * mrSges(6,1) + t8;
t20 = t23 * mrSges(5,3);
t19 = -mrSges(5,1) - t21;
t17 = cos(qJ(3));
t15 = sin(qJ(3));
t13 = t17 ^ 2;
t11 = t15 ^ 2;
t9 = t14 * mrSges(5,2);
t7 = -t16 * pkin(4) - pkin(3);
t5 = t24 * t16;
t4 = -t16 * mrSges(5,1) + t9;
t2 = t24 * t14;
t1 = [t23 * t26 + m(2) + m(3) + m(4); 0; m(3) + m(4) * (t11 + t13) + (t23 * t11 + t13) * t26; m(6) * (t5 * t14 - t2 * t16); (mrSges(4,1) - t3 - t4) * t17 + (t23 * mrSges(6,3) - mrSges(4,2) + t20) * t15 + m(5) * (t23 * t15 * pkin(6) + pkin(3) * t17) + m(6) * (-t7 * t17 + (-t14 * t2 - t16 * t5) * t15); -0.2e1 * pkin(3) * t4 + 0.2e1 * t7 * t3 + Ifges(4,3) + m(6) * (t2 ^ 2 + t5 ^ 2 + t7 ^ 2) + m(5) * (t23 * pkin(6) ^ 2 + pkin(3) ^ 2) + (t5 * t27 + (Ifges(6,2) + Ifges(5,2)) * t16) * t16 + 0.2e1 * pkin(6) * t20 + (t2 * t27 + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t16 + (Ifges(5,1) + Ifges(6,1)) * t14) * t14; t19 * t16 + t8 + t9; ((-mrSges(5,2) - mrSges(6,2)) * t16 + t19 * t14) * t15; t5 * mrSges(6,2) + t21 * t2 + (-mrSges(5,2) * pkin(6) + Ifges(5,6) + Ifges(6,6)) * t16 + (-mrSges(5,1) * pkin(6) - mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t14; Ifges(5,3) + Ifges(6,3) + (0.2e1 * mrSges(6,1) + t25) * pkin(4); 0; -m(6) * t17; m(6) * t7 + t3; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
