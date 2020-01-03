% Return the minimum parameter vector for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% MPV [23x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRPR5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t110 = (pkin(7) ^ 2);
t118 = 2 * pkin(7) * mrSges(6,3) + Ifges(6,2);
t100 = m(6) * t110 + Ifges(5,1) + t118;
t111 = (pkin(4) ^ 2);
t99 = Ifges(5,2) + (t110 + t111) * m(6) + t118;
t120 = t100 - t99;
t119 = (pkin(6) * m(4));
t106 = sin(pkin(9));
t108 = cos(pkin(9));
t117 = t106 * t108;
t103 = t106 ^ 2;
t104 = t108 ^ 2;
t116 = t104 - t103;
t115 = Ifges(5,4) * t117;
t114 = mrSges(4,3) + t119;
t113 = -pkin(7) * m(6) - mrSges(6,3);
t102 = m(6) * pkin(4) + mrSges(5,1);
t112 = -t106 * mrSges(5,2) + t108 * t102;
t109 = cos(pkin(8));
t107 = sin(pkin(8));
t101 = t113 * pkin(4) + Ifges(5,5);
t1 = [Ifges(2,3) + t109 ^ 2 * (pkin(2) ^ 2 * m(4) + Ifges(3,2)) + (0.2e1 * t109 * (t114 * pkin(2) + Ifges(3,4)) + (0.2e1 * t115 + t103 * t100 + t104 * t99 + Ifges(3,1) + Ifges(4,2) + ((2 * mrSges(4,3) + t119) * pkin(6))) * t107) * t107; mrSges(2,1); mrSges(2,2); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t114; mrSges(3,3); m(3) + m(4); t120 * t116 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t115; t116 * Ifges(5,4) + t120 * t117 + Ifges(4,4); -t106 * Ifges(5,6) + t108 * t101 + Ifges(4,5); t108 * Ifges(5,6) + t106 * t101 + Ifges(4,6); (t111 * m(6)) + 0.2e1 * pkin(3) * t112 + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t112; t108 * mrSges(5,2) + t106 * t102 + mrSges(4,2); mrSges(5,3) - t113; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
