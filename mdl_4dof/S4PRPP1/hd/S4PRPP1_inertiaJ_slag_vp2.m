% Calculate joint inertia matrix for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:54
% EndTime: 2019-03-08 18:17:54
% DurationCPUTime: 0.06s
% Computational Cost: add. (21->18), mult. (26->14), div. (0->0), fcn. (0->0), ass. (0->4)
t3 = m(4) + m(5);
t2 = (qJ(3) ^ 2);
t1 = (-pkin(2) - qJ(4));
t4 = [m(2) + m(3) + t3; 0; -2 * pkin(2) * mrSges(4,2) - 2 * t1 * mrSges(5,3) + Ifges(4,1) + Ifges(5,1) + Ifges(3,3) + 2 * (mrSges(4,3) + mrSges(5,2)) * qJ(3) + m(4) * (pkin(2) ^ 2 + t2) + m(5) * (t1 ^ 2 + t2); 0; -m(4) * pkin(2) + m(5) * t1 + mrSges(4,2) - mrSges(5,3); t3; 0; m(5) * qJ(3) + mrSges(5,2); 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t4(1) t4(2) t4(4) t4(7); t4(2) t4(3) t4(5) t4(8); t4(4) t4(5) t4(6) t4(9); t4(7) t4(8) t4(9) t4(10);];
Mq  = res;
