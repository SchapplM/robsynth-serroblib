% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR9
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
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRPR9_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:29
% EndTime: 2019-02-26 21:05:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (39->21), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->31)
t173 = sin(pkin(7));
t177 = cos(pkin(6));
t191 = t173 * t177;
t174 = sin(pkin(6));
t179 = sin(qJ(1));
t190 = t174 * t179;
t181 = cos(qJ(1));
t189 = t174 * t181;
t175 = cos(pkin(12));
t176 = cos(pkin(7));
t188 = t175 * t176;
t172 = sin(pkin(12));
t187 = t179 * t172;
t186 = t179 * t175;
t185 = t181 * t172;
t184 = t181 * t175;
t165 = t177 * t184 - t187;
t183 = -t165 * t176 + t173 * t189;
t167 = -t177 * t186 - t185;
t182 = t167 * t176 + t173 * t190;
t180 = cos(qJ(3));
t178 = sin(qJ(3));
t171 = qJ(4) + pkin(13);
t170 = cos(t171);
t169 = sin(t171);
t168 = -t177 * t187 + t184;
t166 = t177 * t185 + t186;
t164 = -t174 * t175 * t173 + t177 * t176;
t163 = -t167 * t173 + t176 * t190;
t162 = -t165 * t173 - t176 * t189;
t1 = [0, 0, t163, t168 * t178 - t182 * t180, 0 (t168 * t180 + t182 * t178) * t169 - t163 * t170; 0, 0, t162, t166 * t178 + t183 * t180, 0 (t166 * t180 - t183 * t178) * t169 - t162 * t170; 1, 0, t164, -t180 * t191 + (t172 * t178 - t180 * t188) * t174, 0 (t178 * t191 + (t172 * t180 + t178 * t188) * t174) * t169 - t164 * t170;];
Jg_rot  = t1;
