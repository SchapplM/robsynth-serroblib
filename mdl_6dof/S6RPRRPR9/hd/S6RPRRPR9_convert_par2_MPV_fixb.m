% Return the minimum parameter vector for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% MPV [30x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPR9_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t214 = m(4) + m(5);
t237 = (pkin(9) * t214);
t208 = sin(pkin(7));
t226 = mrSges(4,3) + t237;
t236 = t226 * t208;
t212 = cos(pkin(7));
t235 = t226 * t212;
t220 = (pkin(3) ^ 2);
t198 = (t220 * m(5) + Ifges(4,2));
t222 = t198 + (2 * mrSges(4,3) + t237) * pkin(9);
t234 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t206 = sin(pkin(13));
t210 = cos(pkin(13));
t233 = t206 * t210;
t232 = pkin(2) ^ 2 * t214;
t216 = pkin(11) ^ 2;
t192 = m(7) * t216 + Ifges(6,1) + t234;
t219 = (pkin(5) ^ 2);
t196 = m(7) * t219 + Ifges(6,2);
t199 = t206 ^ 2;
t202 = t210 ^ 2;
t230 = t199 * t192 + t202 * t196 + Ifges(5,2);
t229 = pkin(10) * m(5) + mrSges(5,3);
t228 = pkin(11) * m(7) + mrSges(7,3);
t194 = t228 * pkin(5) + Ifges(6,4);
t227 = t194 * t233;
t224 = (2 * pkin(10) * mrSges(5,3)) + 0.2e1 * t227 + t230;
t195 = mrSges(6,2) - t228;
t197 = m(7) * pkin(5) + mrSges(6,1);
t223 = -t206 * t195 + t210 * t197;
t217 = pkin(10) ^ 2;
t213 = cos(pkin(6));
t211 = cos(pkin(12));
t209 = sin(pkin(6));
t207 = sin(pkin(12));
t1 = [Ifges(2,3) + t213 ^ 2 * (t222 * t208 ^ 2 + Ifges(3,3) + t232) + (0.2e1 * (t207 * (-pkin(2) * t235 + Ifges(3,5)) + t211 * (t222 * t212 * t208 + Ifges(3,6))) * t213 + (t211 ^ 2 * (t222 * t212 ^ 2 + Ifges(3,2) + t232) + (0.2e1 * t211 * (pkin(2) * t236 + Ifges(3,4)) + (Ifges(3,1) + t222) * t207) * t207) * t209) * t209; mrSges(2,1); mrSges(2,2); pkin(2) * t214 + mrSges(3,1); mrSges(3,2) - t236; mrSges(3,3) + t235; m(3) + t214; (m(5) * t217) + Ifges(4,1) - t198 + t224; t229 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + ((t217 + t220) * m(5)) + t224; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t229; t192 * t202 + t196 * t199 + Ifges(5,1) - 0.4e1 * t227 - t230; Ifges(5,4) + (t202 - t199) * t194 + (t192 - t196) * t233; Ifges(6,5) * t210 - Ifges(6,6) * t206 + Ifges(5,5); Ifges(6,5) * t206 + Ifges(6,6) * t210 + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t216 + t219) * m(7)) + 0.2e1 * t223 * pkin(4) + t234; mrSges(5,1) + t223; t195 * t210 + t197 * t206 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
