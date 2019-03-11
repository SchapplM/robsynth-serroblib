% Return the minimum parameter vector for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% MPV [31x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRRP8_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t107 = (m(5) + m(6));
t105 = m(4) + t107;
t117 = (pkin(7) * t105);
t116 = (-Ifges(6,2) - Ifges(7,3));
t115 = 2 * pkin(9) * mrSges(6,3) - t116;
t114 = pkin(9) * m(6) + mrSges(6,3);
t113 = -pkin(8) * t107 - mrSges(5,3);
t112 = (mrSges(4,3) - t113);
t111 = (pkin(3) ^ 2);
t110 = (pkin(4) ^ 2);
t109 = (pkin(8) ^ 2);
t108 = pkin(9) ^ 2;
t104 = (t109 + t111);
t1 = [t110 * m(6) + 2 * pkin(8) * mrSges(5,3) + t104 * t107 + Ifges(3,1) + Ifges(4,2) + Ifges(5,2) + Ifges(2,3) + (2 * t112 + t117) * pkin(7); mrSges(2,1); mrSges(2,2); mrSges(3,2) - t112 - t117; mrSges(3,3); m(3) + t105; Ifges(4,1) - Ifges(4,2) + (-t104 + t109) * t107; Ifges(4,4); t113 * pkin(3) + Ifges(4,5); Ifges(4,6); t111 * t107 + Ifges(4,3); pkin(3) * t107 + mrSges(4,1); mrSges(4,2); Ifges(5,1) - Ifges(5,2) + (t108 - t110) * m(6) + t115; t114 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t108 + t110) * m(6) + t115; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t114; Ifges(6,1) + Ifges(7,1) + t116; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
