% Return the minimum parameter vector for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPPPRR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t85 = m(6) + m(7);
t90 = (pkin(7) * t85);
t89 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t88 = pkin(8) * m(7) + mrSges(7,3);
t87 = (pkin(5) ^ 2);
t86 = pkin(8) ^ 2;
t1 = [t87 * m(7) + Ifges(5,1) + Ifges(3,2) + Ifges(6,2) + Ifges(2,3) + Ifges(4,3) + (2 * mrSges(6,3) + t90) * pkin(7); mrSges(2,1); mrSges(2,2); mrSges(3,1); mrSges(3,3); m(3); mrSges(4,1); mrSges(4,2); m(4); mrSges(5,2) - mrSges(6,3) - t90; mrSges(5,3); m(5) + t85; Ifges(6,1) - Ifges(6,2) + (t86 - t87) * m(7) + t89; t88 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t86 + t87) * m(7) + t89; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t88; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
