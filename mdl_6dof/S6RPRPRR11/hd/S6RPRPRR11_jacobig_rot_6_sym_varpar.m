% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRPRR11_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:32
% EndTime: 2019-02-26 20:54:32
% DurationCPUTime: 0.05s
% Computational Cost: add. (39->21), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->31)
t174 = sin(pkin(7));
t178 = cos(pkin(6));
t192 = t174 * t178;
t175 = sin(pkin(6));
t180 = sin(qJ(1));
t191 = t175 * t180;
t182 = cos(qJ(1));
t190 = t175 * t182;
t176 = cos(pkin(12));
t177 = cos(pkin(7));
t189 = t176 * t177;
t173 = sin(pkin(12));
t188 = t180 * t173;
t187 = t180 * t176;
t186 = t182 * t173;
t185 = t182 * t176;
t166 = t178 * t185 - t188;
t184 = -t166 * t177 + t174 * t190;
t168 = -t178 * t187 - t186;
t183 = t168 * t177 + t174 * t191;
t181 = cos(qJ(3));
t179 = sin(qJ(3));
t172 = pkin(13) + qJ(5);
t171 = cos(t172);
t170 = sin(t172);
t169 = -t178 * t188 + t185;
t167 = t178 * t186 + t187;
t165 = -t175 * t176 * t174 + t178 * t177;
t164 = -t168 * t174 + t177 * t191;
t163 = -t166 * t174 - t177 * t190;
t1 = [0, 0, t164, 0, t169 * t179 - t183 * t181 (t169 * t181 + t183 * t179) * t170 - t164 * t171; 0, 0, t163, 0, t167 * t179 + t184 * t181 (t167 * t181 - t184 * t179) * t170 - t163 * t171; 1, 0, t165, 0, -t181 * t192 + (t173 * t179 - t181 * t189) * t175 (t179 * t192 + (t173 * t181 + t179 * t189) * t175) * t170 - t165 * t171;];
Jg_rot  = t1;
