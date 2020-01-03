% Calculate joint inertia matrix for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP4_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:38
% EndTime: 2019-12-31 16:43:38
% DurationCPUTime: 0.16s
% Computational Cost: add. (78->49), mult. (160->60), div. (0->0), fcn. (76->4), ass. (0->15)
t12 = sin(qJ(3));
t13 = cos(qJ(3));
t23 = t12 ^ 2 + t13 ^ 2;
t22 = 0.2e1 * t23;
t21 = -mrSges(4,2) + mrSges(5,3);
t20 = -m(5) * pkin(3) - mrSges(5,1);
t15 = t13 * pkin(3) + t12 * qJ(4);
t11 = cos(pkin(6));
t7 = -t11 * pkin(1) - pkin(2);
t1 = -t15 + t7;
t19 = -0.2e1 * t1;
t10 = sin(pkin(6));
t6 = t10 * pkin(1) + pkin(5);
t18 = t23 * t6 ^ 2;
t2 = [Ifges(2,3) + Ifges(3,3) + (-0.2e1 * t7 * mrSges(4,1) + mrSges(5,1) * t19 + (Ifges(5,3) + Ifges(4,2)) * t13) * t13 + (0.2e1 * t7 * mrSges(4,2) + mrSges(5,3) * t19 + (Ifges(5,1) + Ifges(4,1)) * t12 + 0.2e1 * (Ifges(4,4) - Ifges(5,5)) * t13) * t12 + m(5) * (t1 ^ 2 + t18) + m(4) * (t7 ^ 2 + t18) + (mrSges(5,2) + mrSges(4,3)) * t6 * t22 + (0.2e1 * t11 * mrSges(3,1) - 0.2e1 * t10 * mrSges(3,2) + m(3) * (t10 ^ 2 + t11 ^ 2) * pkin(1)) * pkin(1); 0; m(3) + (m(4) / 0.2e1 + m(5) / 0.2e1) * t22; (qJ(4) * mrSges(5,2) + Ifges(4,6) - Ifges(5,6) + (m(5) * qJ(4) + t21) * t6) * t13 + (-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5) + (-mrSges(4,1) + t20) * t6) * t12; m(5) * t15 + (mrSges(4,1) + mrSges(5,1)) * t13 + t21 * t12; Ifges(5,2) + Ifges(4,3) + 0.2e1 * pkin(3) * mrSges(5,1) + 0.2e1 * qJ(4) * mrSges(5,3) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2); (m(5) * t6 + mrSges(5,2)) * t12; -m(5) * t13; t20; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
