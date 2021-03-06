% Return the minimum parameter vector for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% MPV [32x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPPRR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t148 = (m(3) * pkin(8));
t136 = (m(6) + m(7));
t147 = (pkin(9) * mrSges(6,3));
t140 = (pkin(5) ^ 2);
t146 = (t140 * m(7) + Ifges(6,2));
t145 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t144 = pkin(10) * m(7) + mrSges(7,3);
t143 = pkin(9) * t136 + mrSges(6,3);
t142 = (Ifges(5,1) + Ifges(3,2) + Ifges(4,3) + t146);
t141 = (pkin(4) ^ 2);
t139 = (pkin(9) ^ 2);
t138 = pkin(10) ^ 2;
t135 = sin(pkin(6));
t134 = 2 * t147;
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (t139 * t136 + t134 + (2 * mrSges(3,3) + t148) * pkin(8) + t142) * t135 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t148) * t135; -2 * t147 + Ifges(3,1) + Ifges(4,1) + Ifges(5,2) + (-t139 + t141) * t136 - t142; t143 * pkin(4) + Ifges(3,4) + Ifges(5,4) - Ifges(4,5); Ifges(3,5) + Ifges(4,4) + Ifges(5,6); Ifges(3,6) - Ifges(4,6) + Ifges(5,5); Ifges(3,3) + Ifges(4,2) + Ifges(5,3) + t134 + (t139 + t141) * t136 + t146; mrSges(3,1); mrSges(3,2); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); pkin(4) * t136 + mrSges(5,1); mrSges(5,2) - t143; mrSges(5,3); m(5) + t136; m(7) * t138 + Ifges(6,1) + t145 - t146; t144 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t138 + t140) * m(7) + t145; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t144; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
