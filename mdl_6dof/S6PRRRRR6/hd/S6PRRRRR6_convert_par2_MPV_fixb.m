% Return the minimum parameter vector for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRRRR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_convert_par2_MPV_fixb: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t214 = sin(pkin(8));
t217 = (m(6) + m(7));
t211 = (m(5) + t217);
t227 = pkin(11) * t211 + mrSges(5,3);
t236 = t227 * t214;
t216 = cos(pkin(8));
t235 = t227 * t216;
t234 = (pkin(11) * mrSges(5,3));
t206 = m(4) + t211;
t233 = t206 * pkin(10);
t223 = (pkin(4) ^ 2);
t205 = (t217 * t223 + Ifges(5,2));
t222 = (pkin(5) ^ 2);
t232 = (m(7) * t222 + Ifges(6,2));
t231 = 2 * pkin(13) * mrSges(7,3) + Ifges(7,2);
t230 = 2 * pkin(12) * mrSges(6,3) + t232;
t229 = m(7) * pkin(13) + mrSges(7,3);
t228 = t205 + 2 * t234;
t226 = pkin(12) * t217 + mrSges(6,3);
t221 = (pkin(11) ^ 2);
t225 = t221 * t211 + t228;
t224 = pkin(3) ^ 2;
t220 = pkin(12) ^ 2;
t219 = pkin(13) ^ 2;
t215 = sin(pkin(7));
t210 = t216 ^ 2;
t204 = t210 * t221 + t224;
t203 = mrSges(4,3) + t235;
t1 = [m(2) + m(3) + t206; pkin(2) ^ 2 * t206 + Ifges(3,3) + (t204 * t211 + Ifges(4,2) + t228 * t210 + (0.2e1 * t203 + t233) * pkin(10)) * t215 ^ 2; pkin(2) * t206 + mrSges(3,1); mrSges(3,2) + (-t203 - t233) * t215; -t210 * t205 + Ifges(4,1) - Ifges(4,2) + (-t204 + t221) * t211 + (-0.2e1 * t210 + 0.2e1) * t234 + t205; pkin(3) * t236 + Ifges(4,4); -pkin(3) * t235 + Ifges(4,5); t214 * t216 * t225 + Ifges(4,6); t214 ^ 2 * t225 + t224 * t211 + Ifges(4,3); pkin(3) * t211 + mrSges(4,1); mrSges(4,2) - t236; t217 * t220 + Ifges(5,1) - t205 + t230; pkin(4) * t226 + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t220 + t223) * t217 + t230; pkin(4) * t217 + mrSges(5,1); mrSges(5,2) - t226; m(7) * t219 + Ifges(6,1) + t231 - t232; pkin(5) * t229 + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t219 + t222) * m(7) + t231; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t229; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
