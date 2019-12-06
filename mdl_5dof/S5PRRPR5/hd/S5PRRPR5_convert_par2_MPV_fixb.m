% Return the minimum parameter vector for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% MPV [20x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5PRRPR5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t121 = (pkin(8) ^ 2);
t128 = 2 * pkin(8) * mrSges(6,3) + Ifges(6,2);
t111 = m(6) * t121 + Ifges(5,1) + t128;
t122 = (pkin(4) ^ 2);
t114 = t122 * m(6) + Ifges(5,2);
t129 = t111 - t114;
t119 = sin(pkin(10));
t120 = cos(pkin(10));
t127 = t119 * t120;
t116 = t119 ^ 2;
t117 = t120 ^ 2;
t126 = t117 - t116;
t125 = pkin(8) * m(6) + mrSges(6,3);
t112 = t125 * pkin(4) + Ifges(5,4);
t124 = t112 * t127;
t113 = mrSges(5,2) - t125;
t115 = m(6) * pkin(4) + mrSges(5,1);
t123 = -t119 * t113 + t120 * t115;
t1 = [m(2) + m(3) + m(4); Ifges(3,3) + Ifges(4,2) + t116 * t111 + 0.2e1 * t124 + t117 * t114 + (2 * pkin(7) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * m(4)); m(4) * pkin(2) + mrSges(3,1); -pkin(7) * m(4) + mrSges(3,2) - mrSges(4,3); t129 * t126 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t124; t126 * t112 + t129 * t127 + Ifges(4,4); t120 * Ifges(5,5) - t119 * Ifges(5,6) + Ifges(4,5); t119 * Ifges(5,5) + t120 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t121 + t122) * m(6)) + 0.2e1 * t123 * pkin(3) + t128; mrSges(4,1) + t123; t120 * t113 + t119 * t115 + mrSges(4,2); mrSges(5,3); m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
