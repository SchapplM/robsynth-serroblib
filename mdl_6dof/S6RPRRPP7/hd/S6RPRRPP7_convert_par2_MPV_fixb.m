% Return the minimum parameter vector for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% MPV [28x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPP7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t105 = m(4) + m(5);
t111 = (pkin(7) * t105);
t110 = (-Ifges(5,2) - Ifges(7,2) - Ifges(6,3));
t109 = pkin(8) * m(5) + mrSges(5,3);
t108 = 2 * pkin(8) * mrSges(5,3) - t110;
t107 = (pkin(3) ^ 2);
t106 = pkin(8) ^ 2;
t1 = [t107 * m(5) + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + (2 * mrSges(4,3) + t111) * pkin(7); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t111; mrSges(3,3); m(3) + t105; Ifges(4,1) - Ifges(4,2) + (t106 - t107) * m(5) + t108; t109 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t106 + t107) * m(5) + t108; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t109; Ifges(5,1) + Ifges(6,1) + Ifges(7,1) + t110; Ifges(5,4) - Ifges(6,5) - Ifges(7,4); Ifges(5,5) + Ifges(6,4) - Ifges(7,5); Ifges(5,6) - Ifges(6,6) + Ifges(7,6); Ifges(5,3) + Ifges(6,2) + Ifges(7,3); mrSges(5,1); mrSges(5,2); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
