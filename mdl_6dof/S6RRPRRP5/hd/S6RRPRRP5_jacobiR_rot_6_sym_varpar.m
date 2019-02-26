% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:21
% EndTime: 2019-02-26 21:48:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (154->31), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->37)
t175 = sin(qJ(4));
t179 = cos(qJ(4));
t173 = cos(pkin(6));
t170 = sin(pkin(11));
t172 = cos(pkin(11));
t176 = sin(qJ(2));
t180 = cos(qJ(2));
t183 = t170 * t180 + t172 * t176;
t164 = t183 * t173;
t165 = t170 * t176 - t180 * t172;
t177 = sin(qJ(1));
t181 = cos(qJ(1));
t185 = t181 * t164 - t177 * t165;
t171 = sin(pkin(6));
t189 = t171 * t181;
t150 = t175 * t189 - t179 * t185;
t182 = t165 * t173;
t154 = -t177 * t183 - t181 * t182;
t174 = sin(qJ(5));
t178 = cos(qJ(5));
t195 = t150 * t174 - t154 * t178;
t194 = t150 * t178 + t154 * t174;
t190 = t171 * t177;
t188 = t174 * t179;
t187 = t178 * t179;
t184 = -t177 * t164 - t165 * t181;
t148 = -t175 * t185 - t179 * t189;
t163 = t183 * t171;
t162 = t165 * t171;
t160 = t163 * t179 + t173 * t175;
t159 = -t163 * t175 + t173 * t179;
t157 = t177 * t182 - t181 * t183;
t152 = t175 * t190 + t179 * t184;
t151 = t175 * t184 - t179 * t190;
t147 = t152 * t178 - t157 * t174;
t146 = -t152 * t174 - t157 * t178;
t1 = [t194, t157 * t187 + t174 * t184, 0, -t151 * t178, t146, 0; t147, t154 * t187 + t174 * t185, 0, t148 * t178, t195, 0; 0, -t162 * t187 + t163 * t174, 0, t159 * t178, -t160 * t174 + t162 * t178, 0; -t195, -t157 * t188 + t178 * t184, 0, t151 * t174, -t147, 0; t146, -t154 * t188 + t178 * t185, 0, -t148 * t174, t194, 0; 0, t162 * t188 + t163 * t178, 0, -t159 * t174, -t160 * t178 - t162 * t174, 0; t148, t157 * t175, 0, t152, 0, 0; t151, t154 * t175, 0, -t150, 0, 0; 0, -t162 * t175, 0, t160, 0, 0;];
JR_rot  = t1;
