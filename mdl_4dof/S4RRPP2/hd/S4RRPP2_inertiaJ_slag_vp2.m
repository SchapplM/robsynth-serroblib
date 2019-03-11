% Calculate joint inertia matrix for
% S4RRPP2
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
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:59
% EndTime: 2019-03-08 18:34:00
% DurationCPUTime: 0.10s
% Computational Cost: add. (76->42), mult. (112->43), div. (0->0), fcn. (24->2), ass. (0->16)
t17 = (mrSges(4,3) + mrSges(5,2));
t16 = -2 * mrSges(5,1);
t15 = -mrSges(4,1) - mrSges(5,1);
t13 = Ifges(4,2) + Ifges(3,3) + Ifges(5,3);
t7 = cos(qJ(2));
t5 = -t7 * pkin(1) - pkin(2);
t12 = 2 * t17;
t6 = sin(qJ(2));
t11 = (t7 * mrSges(3,1) - t6 * mrSges(3,2)) * pkin(1);
t9 = qJ(3) ^ 2;
t8 = -pkin(2) - pkin(3);
t4 = pkin(1) * t6 + qJ(3);
t3 = -pkin(3) + t5;
t2 = t4 ^ 2;
t1 = qJ(3) * t4;
t10 = [-0.2e1 * t5 * mrSges(4,1) + t3 * t16 + Ifges(2,3) + t4 * t12 + 0.2e1 * t11 + m(4) * (t5 ^ 2 + t2) + m(5) * (t3 ^ 2 + t2) + m(3) * (t6 ^ 2 + t7 ^ 2) * pkin(1) ^ 2 + t13; t11 + (-t3 - t8) * mrSges(5,1) + (-t5 + pkin(2)) * mrSges(4,1) + m(4) * (-pkin(2) * t5 + t1) + m(5) * (t3 * t8 + t1) + t13 + t17 * (t4 + qJ(3)); 0.2e1 * pkin(2) * mrSges(4,1) + t8 * t16 + qJ(3) * t12 + m(4) * (pkin(2) ^ 2 + t9) + m(5) * (t8 ^ 2 + t9) + t13; m(4) * t5 + m(5) * t3 + t15; -m(4) * pkin(2) + m(5) * t8 + t15; m(4) + m(5); 0; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1) t10(2) t10(4) t10(7); t10(2) t10(3) t10(5) t10(8); t10(4) t10(5) t10(6) t10(9); t10(7) t10(8) t10(9) t10(10);];
Mq  = res;
