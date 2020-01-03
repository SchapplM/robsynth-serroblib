% Calculate joint inertia matrix for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:36
% DurationCPUTime: 0.13s
% Computational Cost: add. (122->48), mult. (211->68), div. (0->0), fcn. (163->4), ass. (0->24)
t18 = sin(pkin(6));
t19 = cos(pkin(6));
t28 = sin(qJ(4));
t29 = cos(qJ(4));
t5 = -t29 * t18 - t28 * t19;
t33 = t5 ^ 2;
t16 = t19 ^ 2;
t20 = -pkin(1) - qJ(3);
t30 = -pkin(5) + t20;
t8 = t30 * t18;
t9 = t30 * t19;
t2 = t28 * t9 + t29 * t8;
t31 = t5 * t2;
t27 = t18 * mrSges(4,1) + t19 * mrSges(4,2);
t26 = t18 ^ 2 + t16;
t7 = -t28 * t18 + t29 * t19;
t25 = t7 ^ 2 + t33;
t24 = -t5 * mrSges(5,1) + t7 * mrSges(5,2);
t23 = m(4) * t26;
t22 = t26 * mrSges(4,3);
t21 = qJ(2) ^ 2;
t10 = t18 * pkin(3) + qJ(2);
t1 = -t28 * t8 + t29 * t9;
t3 = [Ifges(3,1) + Ifges(2,3) + 0.2e1 * t10 * t24 + Ifges(5,2) * t33 - (2 * pkin(1) * mrSges(3,2)) + Ifges(4,1) * t16 + 0.2e1 * mrSges(5,3) * t31 + (-0.2e1 * Ifges(4,4) * t19 + Ifges(4,2) * t18) * t18 - 0.2e1 * t20 * t22 + (-0.2e1 * t1 * mrSges(5,3) + Ifges(5,1) * t7 + 0.2e1 * Ifges(5,4) * t5) * t7 + m(5) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + m(4) * (t26 * t20 ^ 2 + t21) + m(3) * ((pkin(1) ^ 2) + t21) + 0.2e1 * (mrSges(3,3) + t27) * qJ(2); -m(3) * pkin(1) + mrSges(3,2) - t25 * mrSges(5,3) - t22 + m(5) * (t7 * t1 - t31) + t20 * t23; m(5) * t25 + m(3) + t23; m(4) * qJ(2) + m(5) * t10 + t24 + t27; 0; m(4) + m(5); t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t7 + Ifges(5,6) * t5; t7 * mrSges(5,1) + t5 * mrSges(5,2); 0; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
