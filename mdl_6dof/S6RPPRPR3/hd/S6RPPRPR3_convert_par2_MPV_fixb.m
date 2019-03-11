% Return the minimum parameter vector for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% MPV [23x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPPRPR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t103 = sin(pkin(10));
t100 = t103 ^ 2;
t105 = cos(pkin(10));
t101 = t105 ^ 2;
t117 = t101 - t100;
t116 = (pkin(7) * m(5));
t107 = (pkin(8) ^ 2);
t114 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t95 = m(7) * t107 + Ifges(6,1) + t114;
t108 = (pkin(5) ^ 2);
t98 = t108 * m(7) + Ifges(6,2);
t115 = t95 - t98;
t113 = t103 * t105;
t111 = pkin(8) * m(7) + mrSges(7,3);
t96 = t111 * pkin(5) + Ifges(6,4);
t112 = t96 * t113;
t104 = sin(pkin(9));
t106 = cos(pkin(9));
t110 = t106 * mrSges(3,1) - t104 * mrSges(3,2);
t97 = mrSges(6,2) - t111;
t99 = m(7) * pkin(5) + mrSges(6,1);
t109 = -t103 * t97 + t105 * t99;
t1 = [0.2e1 * t112 + t100 * t95 + t101 * t98 + Ifges(4,1) + Ifges(5,2) + Ifges(2,3) + Ifges(3,3) + ((2 * mrSges(5,3) + t116) * pkin(7)) + 0.2e1 * t110 * pkin(1); mrSges(2,1) + t110; t104 * mrSges(3,1) + t106 * mrSges(3,2) + mrSges(2,2); m(3); mrSges(4,2) - mrSges(5,3) - t116; mrSges(4,3); m(4) + m(5); t117 * t115 + Ifges(5,1) - Ifges(5,2) - 0.4e1 * t112; t115 * t113 + t117 * t96 + Ifges(5,4); t105 * Ifges(6,5) - t103 * Ifges(6,6) + Ifges(5,5); t103 * Ifges(6,5) + t105 * Ifges(6,6) + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t107 + t108) * m(7)) + 0.2e1 * t109 * pkin(4) + t114; mrSges(5,1) + t109; t103 * t99 + t105 * t97 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
