% Return the minimum parameter vector for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [26x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPPRPR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t122 = (pkin(8) * mrSges(7,3));
t121 = -pkin(8) * m(7) - mrSges(7,3);
t115 = (pkin(8) ^ 2);
t117 = (pkin(5) ^ 2);
t120 = (Ifges(5,2) + Ifges(7,2) + Ifges(6,3) + (t115 + t117) * m(7));
t112 = sin(pkin(9));
t114 = cos(pkin(9));
t119 = t114 * mrSges(3,1) - t112 * mrSges(3,2);
t109 = 2 * t122;
t118 = 2 * pkin(7) * mrSges(5,3) + t109 + t120;
t116 = pkin(7) ^ 2;
t113 = cos(pkin(10));
t111 = sin(pkin(10));
t1 = [Ifges(2,3) + Ifges(3,3) + t113 ^ 2 * (Ifges(4,2) + (pkin(3) ^ 2 + t116) * m(5) + t118) + (0.2e1 * t113 * Ifges(4,4) + (t116 * m(5) + Ifges(4,1) + t118) * t111) * t111 + 0.2e1 * t119 * pkin(1); mrSges(2,1) + t119; t112 * mrSges(3,1) + t114 * mrSges(3,2) + mrSges(2,2); m(3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); pkin(7) * m(5) + mrSges(4,3) + mrSges(5,3); m(4) + m(5); t117 * m(7) + Ifges(5,1) + Ifges(6,2) - t120 - 2 * t122; Ifges(5,4) + Ifges(6,6); t121 * pkin(5) - Ifges(6,4) + Ifges(5,5); Ifges(5,6) - Ifges(6,5); m(7) * t115 + Ifges(6,1) + Ifges(7,2) + Ifges(5,3) + t109; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) + t121; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
