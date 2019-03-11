% Return the minimum parameter vector for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PPRRRR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_convert_par2_MPV_fixb: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t194 = (m(6) + m(7));
t190 = m(5) + t194;
t205 = (t190 * pkin(10));
t198 = (pkin(5) ^ 2);
t204 = (m(7) * t198 + Ifges(6,2));
t203 = 2 * pkin(12) * mrSges(7,3) + Ifges(7,2);
t202 = 2 * pkin(11) * mrSges(6,3) + t204;
t201 = m(7) * pkin(12) + mrSges(7,3);
t200 = pkin(11) * t194 + mrSges(6,3);
t199 = (pkin(4) ^ 2);
t197 = pkin(11) ^ 2;
t196 = pkin(12) ^ 2;
t193 = sin(pkin(8));
t1 = [m(2); m(3) + m(4) + t190; pkin(3) ^ 2 * t190 + Ifges(4,3) + (t199 * t194 + Ifges(5,2) + (2 * mrSges(5,3) + t205) * pkin(10)) * t193 ^ 2; pkin(3) * t190 + mrSges(4,1); mrSges(4,2) + (-mrSges(5,3) - t205) * t193; Ifges(5,1) - Ifges(5,2) + (t197 - t199) * t194 + t202; pkin(4) * t200 + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t197 + t199) * t194 + t202; pkin(4) * t194 + mrSges(5,1); mrSges(5,2) - t200; m(7) * t196 + Ifges(6,1) + t203 - t204; pkin(5) * t201 + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t196 + t198) * m(7) + t203; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t201; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
