% Calculate joint inertia matrix for
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:41
% DurationCPUTime: 0.12s
% Computational Cost: add. (45->28), mult. (138->43), div. (0->0), fcn. (74->4), ass. (0->15)
t18 = m(5) * pkin(5);
t7 = sin(qJ(4));
t9 = cos(qJ(4));
t2 = -t9 * mrSges(5,1) + t7 * mrSges(5,2);
t17 = m(5) * pkin(3) + mrSges(4,1) - t2;
t5 = t9 ^ 2;
t8 = sin(qJ(3));
t4 = t8 ^ 2;
t10 = cos(qJ(3));
t6 = t10 ^ 2;
t16 = m(3) + m(4) * (t6 + t4);
t15 = t7 ^ 2 + t5;
t14 = mrSges(5,3) * t15;
t12 = -mrSges(5,1) * t7 - mrSges(5,2) * t9;
t1 = [m(2) + m(5) * (t15 * t6 + t4) + t16; m(5) * (-0.1e1 + t15) * t8 * t10; m(5) * (t15 * t4 + t6) + t16; -t17 * t8 + (-mrSges(4,2) + (mrSges(5,3) + t18) * t15) * t10; (t15 * t18 - mrSges(4,2) + t14) * t8 + t17 * t10; Ifges(4,3) + m(5) * (t15 * pkin(5) ^ 2 + pkin(3) ^ 2) + Ifges(5,2) * t5 - 0.2e1 * pkin(3) * t2 + (Ifges(5,1) * t7 + 0.2e1 * Ifges(5,4) * t9) * t7 + 0.2e1 * pkin(5) * t14; t12 * t10; t12 * t8; Ifges(5,5) * t7 + Ifges(5,6) * t9 + t12 * pkin(5); Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
