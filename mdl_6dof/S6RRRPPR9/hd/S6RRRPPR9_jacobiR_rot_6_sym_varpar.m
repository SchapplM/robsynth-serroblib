% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:55
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:54:50
% EndTime: 2019-02-22 11:54:50
% DurationCPUTime: 0.24s
% Computational Cost: add. (156->44), mult. (445->90), div. (0->0), fcn. (628->12), ass. (0->52)
t176 = cos(pkin(6));
t179 = sin(qJ(2));
t184 = cos(qJ(1));
t189 = t184 * t179;
t180 = sin(qJ(1));
t183 = cos(qJ(2));
t191 = t180 * t183;
t168 = t176 * t189 + t191;
t178 = sin(qJ(3));
t182 = cos(qJ(3));
t174 = sin(pkin(6));
t194 = t174 * t184;
t159 = t168 * t182 - t178 * t194;
t188 = t184 * t183;
t192 = t180 * t179;
t167 = -t176 * t188 + t192;
t173 = sin(pkin(11));
t175 = cos(pkin(11));
t146 = t159 * t173 - t167 * t175;
t147 = t159 * t175 + t167 * t173;
t177 = sin(qJ(6));
t181 = cos(qJ(6));
t202 = t146 * t181 - t147 * t177;
t201 = -t146 * t177 - t147 * t181;
t198 = t173 * t182;
t197 = t174 * t178;
t196 = t174 * t182;
t195 = t174 * t183;
t193 = t175 * t182;
t190 = t182 * t183;
t187 = t173 * t181 - t175 * t177;
t186 = t173 * t177 + t175 * t181;
t185 = t168 * t178 + t182 * t194;
t170 = -t176 * t192 + t188;
t169 = t176 * t191 + t189;
t166 = t176 * t178 + t179 * t196;
t165 = t176 * t182 - t179 * t197;
t164 = (t173 * t179 + t175 * t190) * t174;
t163 = (t173 * t190 - t175 * t179) * t174;
t162 = t170 * t182 + t180 * t197;
t161 = -t170 * t178 + t180 * t196;
t157 = t166 * t175 - t173 * t195;
t156 = t166 * t173 + t175 * t195;
t155 = -t169 * t193 + t170 * t173;
t154 = -t169 * t198 - t170 * t175;
t153 = -t167 * t193 + t168 * t173;
t152 = -t167 * t198 - t168 * t175;
t151 = t162 * t175 + t169 * t173;
t150 = t162 * t173 - t169 * t175;
t145 = t150 * t177 + t151 * t181;
t144 = t150 * t181 - t151 * t177;
t1 = [t201, t154 * t177 + t155 * t181, t186 * t161, 0, 0, t144; t145, t152 * t177 + t153 * t181, -t186 * t185, 0, 0, t202; 0, t163 * t177 + t164 * t181, t186 * t165, 0, 0, t156 * t181 - t157 * t177; -t202, t154 * t181 - t155 * t177, t187 * t161, 0, 0, -t145; t144, t152 * t181 - t153 * t177, -t187 * t185, 0, 0, t201; 0, t163 * t181 - t164 * t177, t187 * t165, 0, 0, -t156 * t177 - t157 * t181; t185, t169 * t178, -t162, 0, 0, 0; t161, t167 * t178, -t159, 0, 0, 0; 0, -t178 * t195, -t166, 0, 0, 0;];
JR_rot  = t1;
