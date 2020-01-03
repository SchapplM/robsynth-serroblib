% Calculate joint inertia matrix for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP4_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:46
% EndTime: 2019-12-31 16:27:46
% DurationCPUTime: 0.13s
% Computational Cost: add. (55->43), mult. (125->49), div. (0->0), fcn. (53->2), ass. (0->11)
t7 = sin(qJ(3));
t8 = cos(qJ(3));
t18 = t7 ^ 2 + t8 ^ 2;
t17 = 0.2e1 * t18;
t16 = -mrSges(4,2) + mrSges(5,3);
t15 = -m(5) * pkin(3) - mrSges(5,1);
t10 = t8 * pkin(3) + t7 * qJ(4);
t1 = -pkin(2) - t10;
t14 = -0.2e1 * t1;
t13 = t18 * pkin(5) ^ 2;
t2 = [m(2) + m(3) + (m(4) / 0.2e1 + m(5) / 0.2e1) * t17; 0; Ifges(3,3) + (0.2e1 * pkin(2) * mrSges(4,1) + mrSges(5,1) * t14 + (Ifges(5,3) + Ifges(4,2)) * t8) * t8 + (-0.2e1 * pkin(2) * mrSges(4,2) + mrSges(5,3) * t14 + (Ifges(5,1) + Ifges(4,1)) * t7 + 0.2e1 * (Ifges(4,4) - Ifges(5,5)) * t8) * t7 + m(5) * (t1 ^ 2 + t13) + m(4) * (pkin(2) ^ 2 + t13) + (mrSges(5,2) + mrSges(4,3)) * pkin(5) * t17; m(5) * t10 + (mrSges(4,1) + mrSges(5,1)) * t8 + t16 * t7; (qJ(4) * mrSges(5,2) + Ifges(4,6) - Ifges(5,6) + (m(5) * qJ(4) + t16) * pkin(5)) * t8 + (-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5) + (-mrSges(4,1) + t15) * pkin(5)) * t7; Ifges(5,2) + Ifges(4,3) + 0.2e1 * pkin(3) * mrSges(5,1) + 0.2e1 * qJ(4) * mrSges(5,3) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2); -m(5) * t8; (m(5) * pkin(5) + mrSges(5,2)) * t7; t15; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
