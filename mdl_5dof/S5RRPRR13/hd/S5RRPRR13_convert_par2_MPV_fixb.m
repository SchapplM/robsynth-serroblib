% Return the minimum parameter vector for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% MPV [28x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPRR13_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t124 = -pkin(8) * m(6) - mrSges(6,3);
t104 = (mrSges(5,3) - t124);
t111 = (m(5) + m(6));
t123 = -pkin(7) * t111 - t104;
t113 = (pkin(8) ^ 2);
t115 = (pkin(4) ^ 2);
t121 = (Ifges(5,2) + (t113 + t115) * m(6));
t109 = sin(pkin(9));
t110 = cos(pkin(9));
t120 = t109 * t110;
t119 = Ifges(4,4) * t120;
t116 = (pkin(3) ^ 2);
t118 = (t116 * t111 + Ifges(3,2) + Ifges(4,3));
t117 = 2 * pkin(8) * mrSges(6,3) + 2 * pkin(7) * t104 + Ifges(6,2) + t121;
t114 = pkin(7) ^ 2;
t107 = t110 ^ 2;
t106 = t109 ^ 2;
t101 = t123 * pkin(3) + Ifges(4,5);
t100 = t114 * t111 + Ifges(4,1) + t117;
t99 = Ifges(4,2) + (t114 + t116) * t111 + t117;
t1 = [Ifges(2,3) + 2 * pkin(6) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(6) ^ 2) * m(3) + t118; m(3) * pkin(1) + mrSges(2,1); -pkin(6) * m(3) + mrSges(2,2) - mrSges(3,3); t107 * t100 + t106 * t99 + Ifges(3,1) - t118 - 0.2e1 * t119; t109 * Ifges(4,6) - t110 * t101 + Ifges(3,4); Ifges(3,5) + (t107 - t106) * Ifges(4,4) + (t100 - t99) * t120; -t110 * Ifges(4,6) - t109 * t101 + Ifges(3,6); t106 * t100 + t107 * t99 + Ifges(3,3) + 0.2e1 * t119; mrSges(3,1); mrSges(3,2); pkin(3) * t111 + mrSges(4,1); mrSges(4,2); mrSges(4,3) - t123; m(4) + t111; m(6) * t113 + Ifges(5,1) - t121; Ifges(5,4); t124 * pkin(4) + Ifges(5,5); Ifges(5,6); t115 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
