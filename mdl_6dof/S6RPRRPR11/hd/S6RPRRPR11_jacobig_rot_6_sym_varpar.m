% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRPR11_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:44
% EndTime: 2019-02-26 21:06:44
% DurationCPUTime: 0.07s
% Computational Cost: add. (33->20), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->30)
t175 = sin(pkin(7));
t179 = cos(pkin(6));
t195 = t175 * t179;
t176 = sin(pkin(6));
t182 = sin(qJ(1));
t194 = t176 * t182;
t185 = cos(qJ(1));
t193 = t176 * t185;
t177 = cos(pkin(12));
t178 = cos(pkin(7));
t192 = t177 * t178;
t174 = sin(pkin(12));
t191 = t182 * t174;
t190 = t182 * t177;
t189 = t185 * t174;
t188 = t185 * t177;
t170 = t179 * t188 - t191;
t187 = -t170 * t178 + t175 * t193;
t172 = -t179 * t190 - t189;
t186 = t172 * t178 + t175 * t194;
t184 = cos(qJ(3));
t183 = cos(qJ(4));
t181 = sin(qJ(3));
t180 = sin(qJ(4));
t173 = -t179 * t191 + t188;
t171 = t179 * t189 + t190;
t169 = -t176 * t177 * t175 + t179 * t178;
t168 = -t172 * t175 + t178 * t194;
t167 = -t170 * t175 - t178 * t193;
t1 = [0, 0, t168, t173 * t181 - t186 * t184, 0 (t173 * t184 + t186 * t181) * t180 - t168 * t183; 0, 0, t167, t171 * t181 + t187 * t184, 0 (t171 * t184 - t187 * t181) * t180 - t167 * t183; 1, 0, t169, -t184 * t195 + (t174 * t181 - t184 * t192) * t176, 0 (t181 * t195 + (t174 * t184 + t181 * t192) * t176) * t180 - t169 * t183;];
Jg_rot  = t1;
