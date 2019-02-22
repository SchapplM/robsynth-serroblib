% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:14
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:14:26
% EndTime: 2019-02-22 11:14:26
% DurationCPUTime: 0.16s
% Computational Cost: add. (157->33), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->37)
t170 = sin(pkin(11));
t172 = cos(pkin(11));
t176 = sin(qJ(2));
t180 = cos(qJ(2));
t165 = t176 * t170 - t180 * t172;
t173 = cos(pkin(6));
t163 = t165 * t173;
t177 = sin(qJ(1));
t181 = cos(qJ(1));
t182 = t180 * t170 + t176 * t172;
t154 = -t181 * t163 - t177 * t182;
t175 = sin(qJ(5));
t179 = cos(qJ(5));
t171 = sin(pkin(6));
t189 = t171 * t181;
t150 = t154 * t175 + t179 * t189;
t174 = sin(qJ(6));
t178 = cos(qJ(6));
t164 = t182 * t173;
t184 = t181 * t164 - t177 * t165;
t194 = t150 * t174 + t184 * t178;
t193 = t150 * t178 - t184 * t174;
t190 = t171 * t177;
t188 = t174 * t175;
t187 = t175 * t178;
t149 = -t154 * t179 + t175 * t189;
t183 = -t177 * t164 - t181 * t165;
t162 = t182 * t171;
t161 = t165 * t171;
t160 = t161 * t175 + t173 * t179;
t159 = t161 * t179 - t173 * t175;
t157 = t177 * t163 - t181 * t182;
t148 = -t157 * t175 + t179 * t190;
t147 = t157 * t179 + t175 * t190;
t146 = t148 * t178 + t174 * t183;
t145 = -t148 * t174 + t178 * t183;
t1 = [t193, t157 * t174 + t183 * t187, 0, 0, -t147 * t178, t145; t146, t154 * t174 + t184 * t187, 0, 0, t149 * t178, t194; 0, -t161 * t174 + t162 * t187, 0, 0, t159 * t178, -t160 * t174 + t162 * t178; -t194, t157 * t178 - t183 * t188, 0, 0, t147 * t174, -t146; t145, t154 * t178 - t184 * t188, 0, 0, -t149 * t174, t193; 0, -t161 * t178 - t162 * t188, 0, 0, -t159 * t174, -t160 * t178 - t162 * t174; t149, -t183 * t179, 0, 0, t148, 0; t147, -t184 * t179, 0, 0, -t150, 0; 0, -t162 * t179, 0, 0, t160, 0;];
JR_rot  = t1;
